#include <memory>
#include <string>

#include <gz/common/Console.hh>
#include <gz/plugin/Register.hh>
#include <gz/sim/Entity.hh>
#include <gz/sim/EntityComponentManager.hh>
#include <gz/sim/EventManager.hh>
#include <gz/sim/Link.hh>
#include <gz/sim/Model.hh>
#include <gz/sim/System.hh>
#include <sdf/Element.hh>

#include "usv_hydro/BuoyancyModel.hh"
#include "usv_hydro/DragModel.hh"
#include "usv_hydro/EnvironmentModel.hh"
#include "usv_hydro/HydroIntegrator.hh"
#include "usv_hydro/HydroTypes.hh"

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

    this->config = this->LoadConfig(_sdf);

    std::string configError;
    if (!this->config.IsValid(&configError))
    {
      gzerr << "[usv_hydro] Invalid plugin config: " << configError << "\n";
      return;
    }

    this->environment.SetFromConfig(this->config);
    this->buoyancy.SetParameters(
        this->config.fluidDensity, this->config.gravity);
    this->drag.SetParameters(
        this->config.fluidDensity, this->config.cd, this->config.linearDrag);

    auto linkEntity = this->model.LinkByName(_ecm, this->config.linkName);
    if (linkEntity == gz::sim::kNullEntity)
      linkEntity = this->model.CanonicalLink(_ecm);

    if (linkEntity == gz::sim::kNullEntity)
    {
      gzerr << "[usv_hydro] Could not find link [" << this->config.linkName
            << "] or canonical link.\n";
      return;
    }

    this->link = gz::sim::Link(linkEntity);
    this->link.EnableVelocityChecks(_ecm, true);
    this->link.EnableBoundingBoxChecks(_ecm, true);

    this->cellsReady = false;
    this->cellGrid = HydroCellGrid{};
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

    if (this->cellGrid.offsets.empty())
      return;

    const HydroKinematics kinematics{
        pose->Pos(),
        pose->Rot(),
        *linearVel,
        *angularVel};

    for (const auto &offsetLocal : this->cellGrid.offsets)
    {
      const auto cellForces = this->integrator.ComputeCellForces(
          offsetLocal,
          this->cellGrid,
          kinematics,
          this->environment,
          this->buoyancy,
          this->drag);

      if (!cellForces.submerged)
        continue;

      this->link.AddWorldForce(_ecm, cellForces.buoyancyWorld, offsetLocal);
      this->link.AddWorldForce(_ecm, cellForces.dragWorld, offsetLocal);
    }
  }

  private: HydroConfig LoadConfig(
      const std::shared_ptr<const sdf::Element> &_sdf) const
  {
    HydroConfig result;
    if (!_sdf)
      return result;

    if (_sdf->HasElement("link_name"))
      result.linkName = _sdf->Get<std::string>("link_name");
    if (_sdf->HasElement("fluid_density"))
      result.fluidDensity = _sdf->Get<double>("fluid_density");
    if (_sdf->HasElement("water_level"))
      result.waterLevel = _sdf->Get<double>("water_level");
    if (_sdf->HasElement("gravity"))
      result.gravity = _sdf->Get<double>("gravity");

    if (_sdf->HasElement("cells_x"))
      result.cellsX = _sdf->Get<int>("cells_x");
    if (_sdf->HasElement("cells_y"))
      result.cellsY = _sdf->Get<int>("cells_y");
    if (_sdf->HasElement("cells_z"))
      result.cellsZ = _sdf->Get<int>("cells_z");

    if (_sdf->HasElement("cd_x"))
      result.cd.X(_sdf->Get<double>("cd_x"));
    if (_sdf->HasElement("cd_y"))
      result.cd.Y(_sdf->Get<double>("cd_y"));
    if (_sdf->HasElement("cd_z"))
      result.cd.Z(_sdf->Get<double>("cd_z"));

    if (_sdf->HasElement("linear_drag_x"))
      result.linearDrag.X(_sdf->Get<double>("linear_drag_x"));
    if (_sdf->HasElement("linear_drag_y"))
      result.linearDrag.Y(_sdf->Get<double>("linear_drag_y"));
    if (_sdf->HasElement("linear_drag_z"))
      result.linearDrag.Z(_sdf->Get<double>("linear_drag_z"));

    if (_sdf->HasElement("current_velocity"))
      result.currentVelocity = _sdf->Get<gz::math::Vector3d>("current_velocity");

    return result;
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

    const int count = this->config.cellsX * this->config.cellsY * this->config.cellsZ;
    if (count <= 0)
      return false;

    this->cellGrid = HydroCellGrid{};
    this->cellGrid.offsets.reserve(static_cast<std::size_t>(count));

    for (int ix = 0; ix < this->config.cellsX; ++ix)
    {
      for (int iy = 0; iy < this->config.cellsY; ++iy)
      {
        for (int iz = 0; iz < this->config.cellsZ; ++iz)
        {
          const double x = min.X() + (ix + 0.5) * size.X() / this->config.cellsX;
          const double y = min.Y() + (iy + 0.5) * size.Y() / this->config.cellsY;
          const double z = min.Z() + (iz + 0.5) * size.Z() / this->config.cellsZ;
          this->cellGrid.offsets.emplace_back(x, y, z);
        }
      }
    }

    this->cellGrid.cellVolume = (size.X() * size.Y() * size.Z()) / count;
    this->cellGrid.cellHeightApprox = size.Z() / this->config.cellsZ;
    this->cellGrid.dragAreaCell.Set(
        (size.Y() * size.Z()) / count,
        (size.X() * size.Z()) / count,
        (size.X() * size.Y()) / count);

    this->cellsReady = true;
    gzmsg << "[usv_hydro] Hydro cells initialized for link ["
          << this->config.linkName << "] count=" << count << "\n";
    return true;
  }

  private: bool configured{false};
  private: bool cellsReady{false};

  private: HydroConfig config;
  private: HydroCellGrid cellGrid;
  private: EnvironmentModel environment;
  private: BuoyancyModel buoyancy;
  private: DragModel drag;
  private: HydroIntegrator integrator;

  private: gz::sim::Model model;
  private: gz::sim::Link link;
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
