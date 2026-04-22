#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <gz/common/Console.hh>
#include <gz/math/Pose3.hh>
#include <gz/math/Quaternion.hh>
#include <gz/math/Vector3.hh>
#include <gz/plugin/Register.hh>
#include <gz/sim/Entity.hh>
#include <gz/sim/EntityComponentManager.hh>
#include <gz/sim/EventManager.hh>
#include <gz/sim/Link.hh>
#include <gz/sim/Model.hh>
#include <gz/sim/System.hh>
#include <sdf/Element.hh>

namespace usv_hydro
{
namespace systems
{

class HydroCellSystem final:
  public gz::sim::System,
  public gz::sim::ISystemConfigure,
  public gz::sim::ISystemPreUpdate
{
  public: void Configure(
      const gz::sim::Entity &_entity,
      const std::shared_ptr<const sdf::Element> &_sdf,
      gz::sim::EntityComponentManager &_ecm,
      gz::sim::EventManager & /*_eventMgr*/) override
  {
    this->model = gz::sim::Model(_entity);
    if (!this->model.Valid(_ecm))
    {
      gzerr << "[usv_hydro] Plugin must be attached to a model.\n";
      return;
    }

    if (_sdf->HasElement("link_name"))
      this->linkName = _sdf->Get<std::string>("link_name");
    if (_sdf->HasElement("fluid_density"))
      this->fluidDensity = _sdf->Get<double>("fluid_density");
    if (_sdf->HasElement("water_level"))
      this->waterLevel = _sdf->Get<double>("water_level");
    if (_sdf->HasElement("gravity"))
      this->gravity = _sdf->Get<double>("gravity");

    if (_sdf->HasElement("cells_x"))
      this->cellsX = std::max(1, _sdf->Get<int>("cells_x"));
    if (_sdf->HasElement("cells_y"))
      this->cellsY = std::max(1, _sdf->Get<int>("cells_y"));
    if (_sdf->HasElement("cells_z"))
      this->cellsZ = std::max(1, _sdf->Get<int>("cells_z"));

    if (_sdf->HasElement("cd_x"))
      this->cd.X(_sdf->Get<double>("cd_x"));
    if (_sdf->HasElement("cd_y"))
      this->cd.Y(_sdf->Get<double>("cd_y"));
    if (_sdf->HasElement("cd_z"))
      this->cd.Z(_sdf->Get<double>("cd_z"));

    if (_sdf->HasElement("linear_drag_x"))
      this->linearDrag.X(_sdf->Get<double>("linear_drag_x"));
    if (_sdf->HasElement("linear_drag_y"))
      this->linearDrag.Y(_sdf->Get<double>("linear_drag_y"));
    if (_sdf->HasElement("linear_drag_z"))
      this->linearDrag.Z(_sdf->Get<double>("linear_drag_z"));

    if (_sdf->HasElement("current_velocity"))
      this->currentVelocity = _sdf->Get<gz::math::Vector3d>("current_velocity");

    auto linkEntity = this->model.LinkByName(_ecm, this->linkName);
    if (linkEntity == gz::sim::kNullEntity)
    {
      linkEntity = this->model.CanonicalLink(_ecm);
    }

    if (linkEntity == gz::sim::kNullEntity)
    {
      gzerr << "[usv_hydro] Could not find link [" << this->linkName
            << "] or canonical link.\n";
      return;
    }

    this->link = gz::sim::Link(linkEntity);
    this->link.EnableVelocityChecks(_ecm, true);
    this->link.EnableBoundingBoxChecks(_ecm, true);
    this->configured = true;
  }

  public: void PreUpdate(
      const gz::sim::UpdateInfo &_info,
      gz::sim::EntityComponentManager &_ecm) override
  {
    if (!this->configured || _info.paused)
      return;

    if (!this->cellsReady && !this->InitializeCells(_ecm))
      return;

    const auto pose = this->link.WorldPose(_ecm);
    const auto linearVel = this->link.WorldLinearVelocity(_ecm);
    const auto angularVel = this->link.WorldAngularVelocity(_ecm);
    if (!pose || !linearVel || !angularVel)
      return;

    const auto &rot = pose->Rot();
    const auto nCells = static_cast<int>(this->cellOffsets.size());
    if (nCells <= 0)
      return;

    for (const auto &offsetLocal : this->cellOffsets)
    {
      const auto offsetWorld = rot.RotateVector(offsetLocal);
      const auto cellPosWorld = pose->Pos() + offsetWorld;

      const double depth = this->waterLevel - cellPosWorld.Z();
      if (depth <= 0.0)
        continue;

      const double submergence =
          std::clamp(depth / this->cellHeightApprox, 0.0, 1.0);
      const double displacedVolume = this->cellVolume * submergence;

      const gz::math::Vector3d buoyancyWorld(
          0.0, 0.0, this->fluidDensity * this->gravity * displacedVolume);
      this->link.AddWorldForce(_ecm, buoyancyWorld, offsetLocal);

      const auto cellVelWorld =
          *linearVel + angularVel->Cross(offsetWorld) - this->currentVelocity;
      const auto cellVelBody = rot.RotateVectorReverse(cellVelWorld);

      gz::math::Vector3d dragBody;
      dragBody.X(
          -this->linearDrag.X() * cellVelBody.X()
          - 0.5 * this->fluidDensity * this->cd.X() * this->dragAreaCell.X()
          * std::abs(cellVelBody.X()) * cellVelBody.X());
      dragBody.Y(
          -this->linearDrag.Y() * cellVelBody.Y()
          - 0.5 * this->fluidDensity * this->cd.Y() * this->dragAreaCell.Y()
          * std::abs(cellVelBody.Y()) * cellVelBody.Y());
      dragBody.Z(
          -this->linearDrag.Z() * cellVelBody.Z()
          - 0.5 * this->fluidDensity * this->cd.Z() * this->dragAreaCell.Z()
          * std::abs(cellVelBody.Z()) * cellVelBody.Z());

      const auto dragWorld = rot.RotateVector(dragBody * submergence);
      this->link.AddWorldForce(_ecm, dragWorld, offsetLocal);
    }
  }

  private: bool InitializeCells(gz::sim::EntityComponentManager &_ecm)
  {
    const auto aabb = this->link.AxisAlignedBox(_ecm);
    if (!aabb)
      return false;

    const auto min = aabb->Min();
    const auto max = aabb->Max();
    const auto size = max - min;

    if (size.X() <= 1e-6 || size.Y() <= 1e-6 || size.Z() <= 1e-6)
    {
      gzerr << "[usv_hydro] Invalid hull AABB dimensions.\n";
      return false;
    }

    const int count = this->cellsX * this->cellsY * this->cellsZ;
    this->cellOffsets.reserve(count);

    for (int ix = 0; ix < this->cellsX; ++ix)
    {
      for (int iy = 0; iy < this->cellsY; ++iy)
      {
        for (int iz = 0; iz < this->cellsZ; ++iz)
        {
          const double x = min.X() + (ix + 0.5) * size.X() / this->cellsX;
          const double y = min.Y() + (iy + 0.5) * size.Y() / this->cellsY;
          const double z = min.Z() + (iz + 0.5) * size.Z() / this->cellsZ;
          this->cellOffsets.emplace_back(x, y, z);
        }
      }
    }

    this->cellVolume = (size.X() * size.Y() * size.Z()) / count;
    this->cellHeightApprox = size.Z() / this->cellsZ;
    this->dragAreaCell.Set(
        (size.Y() * size.Z()) / count,
        (size.X() * size.Z()) / count,
        (size.X() * size.Y()) / count);

    this->cellsReady = true;
    gzmsg << "[usv_hydro] Hydro cells initialized for link [" << this->linkName
          << "] count=" << count << "\n";
    return true;
  }

  private: bool configured{false};
  private: bool cellsReady{false};
  private: std::string linkName{"hull_link"};
  private: gz::sim::Model model;
  private: gz::sim::Link link;

  private: double fluidDensity{1000.0};
  private: double waterLevel{0.0};
  private: double gravity{9.81};

  private: int cellsX{12};
  private: int cellsY{6};
  private: int cellsZ{4};
  private: std::vector<gz::math::Vector3d> cellOffsets;

  private: double cellVolume{0.0};
  private: double cellHeightApprox{0.1};
  private: gz::math::Vector3d dragAreaCell{0.0, 0.0, 0.0};

  private: gz::math::Vector3d cd{0.9, 1.2, 2.0};
  private: gz::math::Vector3d linearDrag{20.0, 40.0, 70.0};
  private: gz::math::Vector3d currentVelocity{0.0, 0.0, 0.0};
};

}  // namespace systems
}  // namespace usv_hydro

GZ_ADD_PLUGIN(
    usv_hydro::systems::HydroCellSystem,
    gz::sim::System,
    gz::sim::ISystemConfigure,
    gz::sim::ISystemPreUpdate)

GZ_ADD_PLUGIN_ALIAS(
    usv_hydro::systems::HydroCellSystem,
    "usv_hydro::systems::HydroCellSystem")
