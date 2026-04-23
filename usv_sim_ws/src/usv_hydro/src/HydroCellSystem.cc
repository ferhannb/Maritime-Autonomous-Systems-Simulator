#include <algorithm>
#include <cmath>
#include <cstdlib>
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
  private: struct HydrostaticAggregate
  {
    gz::math::Vector3d forceWorld{0.0, 0.0, 0.0};
    gz::math::Vector3d momentWorld{0.0, 0.0, 0.0};
  };

  private: struct LinearHydrostaticStiffness
  {
    double k33{0.0};
    double k34{0.0};
    double k35{0.0};
    double k43{0.0};
    double k44{0.0};
    double k45{0.0};
    double k53{0.0};
    double k54{0.0};
    double k55{0.0};
  };

  private: struct LinearHydrostaticWrench
  {
    gz::math::Vector3d forceWorld{0.0, 0.0, 0.0};
    gz::math::Vector3d torqueWorld{0.0, 0.0, 0.0};
  };

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
    this->useHydrostaticStiffnessMatrix = this->config.useHydrostaticStiffnessMatrix;
    this->hydrostaticStiffnessScale = this->config.hydrostaticStiffnessScale;
    this->stiffnessHeaveStep = this->config.stiffnessHeaveStep;
    this->stiffnessAngleStep = this->config.stiffnessAngleStep;

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
    this->hydrostaticStiffnessReady = false;
    this->hydrostaticStiffness = LinearHydrostaticStiffness{};
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

    if (this->useHydrostaticStiffnessMatrix &&
        !this->hydrostaticStiffnessReady)
    {
      if (!this->EstimateHydrostaticStiffness(*pose))
      {
        gzerr << "[usv_hydro] Failed to estimate hydrostatic stiffness matrix; "
              << "disabling matrix runtime contribution.\n";
        this->useHydrostaticStiffnessMatrix = false;
      }
      else
      {
        this->hydrostaticStiffnessReady = true;
        gzmsg << "[usv_hydro] Runtime hydrostatic stiffness enabled "
              << "(scale=" << this->hydrostaticStiffnessScale << ") "
              << "K33=" << this->hydrostaticStiffness.k33
              << " K44=" << this->hydrostaticStiffness.k44
              << " K55=" << this->hydrostaticStiffness.k55 << "\n";
      }
    }

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

    if (this->useHydrostaticStiffnessMatrix &&
        this->hydrostaticStiffnessReady &&
        this->hydrostaticStiffnessScale > 0.0)
    {
      const auto restoringWrench = this->ComputeLinearHydrostaticWrench(*pose);
      this->link.AddWorldWrench(
          _ecm, restoringWrench.forceWorld, restoringWrench.torqueWorld);
    }
  }

  private: HydrostaticAggregate ComputeHydrostaticAggregate(
      const gz::math::Vector3d &_positionWorld,
      const gz::math::Quaterniond &_rotationWorld) const
  {
    HydrostaticAggregate aggregate;

    for (const auto &offsetLocal : this->cellGrid.offsets)
    {
      const auto offsetWorld = _rotationWorld.RotateVector(offsetLocal);
      const auto cellPositionWorld = _positionWorld + offsetWorld;
      const double depth = this->environment.DepthAt(cellPositionWorld.Z());
      const double submergence =
          this->buoyancy.ComputeSubmergence(depth, this->cellGrid.cellHeightApprox);
      if (submergence <= 0.0)
        continue;

      const auto buoyancyForce =
          this->buoyancy.ComputeForceWorld(this->cellGrid.cellVolume, submergence);
      aggregate.forceWorld += buoyancyForce;
      aggregate.momentWorld += offsetWorld.Cross(buoyancyForce);
    }

    return aggregate;
  }

  private: bool EstimateHydrostaticStiffness(const gz::math::Pose3d &_referencePose)
  {
    if (this->cellGrid.offsets.empty())
      return false;

    this->referencePositionWorld = _referencePose.Pos();
    this->referenceRotationWorld = _referencePose.Rot();

    const auto aggregateAt = [&](const double _dz,
                                 const double _roll,
                                 const double _pitch)
    {
      const auto positionWorld =
          this->referencePositionWorld +
          this->referenceRotationWorld.RotateVector(
              gz::math::Vector3d(0.0, 0.0, _dz));
      const auto rotationWorld =
          this->referenceRotationWorld * gz::math::Quaterniond(_roll, _pitch, 0.0);
      return this->ComputeHydrostaticAggregate(positionWorld, rotationWorld);
    };

    const double dz = this->stiffnessHeaveStep;
    const double dTheta = this->stiffnessAngleStep;

    const auto fZPlus = aggregateAt(dz, 0.0, 0.0).forceWorld.Z();
    const auto fZMinus = aggregateAt(-dz, 0.0, 0.0).forceWorld.Z();
    const auto fZRollPlus = aggregateAt(0.0, dTheta, 0.0).forceWorld.Z();
    const auto fZRollMinus = aggregateAt(0.0, -dTheta, 0.0).forceWorld.Z();
    const auto fZPitchPlus = aggregateAt(0.0, 0.0, dTheta).forceWorld.Z();
    const auto fZPitchMinus = aggregateAt(0.0, 0.0, -dTheta).forceWorld.Z();

    const auto mXPlus = aggregateAt(0.0, dTheta, 0.0).momentWorld.X();
    const auto mXMinus = aggregateAt(0.0, -dTheta, 0.0).momentWorld.X();
    const auto mXHeavePlus = aggregateAt(dz, 0.0, 0.0).momentWorld.X();
    const auto mXHeaveMinus = aggregateAt(-dz, 0.0, 0.0).momentWorld.X();
    const auto mXPitchPlus = aggregateAt(0.0, 0.0, dTheta).momentWorld.X();
    const auto mXPitchMinus = aggregateAt(0.0, 0.0, -dTheta).momentWorld.X();

    const auto mYPlus = aggregateAt(0.0, 0.0, dTheta).momentWorld.Y();
    const auto mYMinus = aggregateAt(0.0, 0.0, -dTheta).momentWorld.Y();
    const auto mYHeavePlus = aggregateAt(dz, 0.0, 0.0).momentWorld.Y();
    const auto mYHeaveMinus = aggregateAt(-dz, 0.0, 0.0).momentWorld.Y();
    const auto mYRollPlus = aggregateAt(0.0, dTheta, 0.0).momentWorld.Y();
    const auto mYRollMinus = aggregateAt(0.0, -dTheta, 0.0).momentWorld.Y();

    this->hydrostaticStiffness.k33 = -(fZPlus - fZMinus) / (2.0 * dz);
    this->hydrostaticStiffness.k34 = -(fZRollPlus - fZRollMinus) / (2.0 * dTheta);
    this->hydrostaticStiffness.k35 = -(fZPitchPlus - fZPitchMinus) / (2.0 * dTheta);
    this->hydrostaticStiffness.k43 = -(mXHeavePlus - mXHeaveMinus) / (2.0 * dz);
    this->hydrostaticStiffness.k44 = -(mXPlus - mXMinus) / (2.0 * dTheta);
    this->hydrostaticStiffness.k45 = -(mXPitchPlus - mXPitchMinus) / (2.0 * dTheta);
    this->hydrostaticStiffness.k53 = -(mYHeavePlus - mYHeaveMinus) / (2.0 * dz);
    this->hydrostaticStiffness.k54 = -(mYRollPlus - mYRollMinus) / (2.0 * dTheta);
    this->hydrostaticStiffness.k55 = -(mYPlus - mYMinus) / (2.0 * dTheta);

    return std::isfinite(this->hydrostaticStiffness.k33) &&
           std::isfinite(this->hydrostaticStiffness.k34) &&
           std::isfinite(this->hydrostaticStiffness.k35) &&
           std::isfinite(this->hydrostaticStiffness.k43) &&
           std::isfinite(this->hydrostaticStiffness.k44) &&
           std::isfinite(this->hydrostaticStiffness.k45) &&
           std::isfinite(this->hydrostaticStiffness.k53) &&
           std::isfinite(this->hydrostaticStiffness.k54) &&
           std::isfinite(this->hydrostaticStiffness.k55);
  }

  private: LinearHydrostaticWrench ComputeLinearHydrostaticWrench(
      const gz::math::Pose3d &_pose) const
  {
    LinearHydrostaticWrench wrench;

    const auto deltaPositionWorld = _pose.Pos() - this->referencePositionWorld;
    const auto deltaPositionRef =
        this->referenceRotationWorld.RotateVectorReverse(deltaPositionWorld);
    const auto rotationError =
        this->referenceRotationWorld.Inverse() * _pose.Rot();
    const auto eulerError = rotationError.Euler();

    const double dz = deltaPositionRef.Z();
    const double roll = eulerError.X();
    const double pitch = eulerError.Y();
    const double scale = this->hydrostaticStiffnessScale;

    const gz::math::Vector3d forceRef(
        0.0,
        0.0,
        -scale *
        (this->hydrostaticStiffness.k33 * dz +
         this->hydrostaticStiffness.k34 * roll +
         this->hydrostaticStiffness.k35 * pitch));

    const gz::math::Vector3d torqueRef(
        -scale *
        (this->hydrostaticStiffness.k43 * dz +
         this->hydrostaticStiffness.k44 * roll +
         this->hydrostaticStiffness.k45 * pitch),
        -scale *
        (this->hydrostaticStiffness.k53 * dz +
         this->hydrostaticStiffness.k54 * roll +
         this->hydrostaticStiffness.k55 * pitch),
        0.0);

    wrench.forceWorld = this->referenceRotationWorld.RotateVector(forceRef);
    wrench.torqueWorld = this->referenceRotationWorld.RotateVector(torqueRef);
    return wrench;
  }

  private: HydroConfig LoadConfig(
      const std::shared_ptr<const sdf::Element> &_sdf) const
  {
    HydroConfig result;
    if (!_sdf)
      return result;

    if (_sdf->HasElement("link_name"))
      result.linkName = _sdf->Get<std::string>("link_name");

    std::string configFile;
    if (_sdf->HasElement("config_file"))
      configFile = _sdf->Get<std::string>("config_file");
    if (_sdf->HasElement("profile"))
      result.profile = _sdf->Get<std::string>("profile");

    if (const char *envConfig = std::getenv("USV_HYDRO_CONFIG_FILE"))
    {
      if (envConfig[0] != '\0')
      {
        configFile = envConfig;
        gzmsg << "[usv_hydro] Using config_file from USV_HYDRO_CONFIG_FILE: "
              << configFile << "\n";
      }
    }
    if (const char *envProfile = std::getenv("USV_HYDRO_PROFILE"))
    {
      if (envProfile[0] != '\0')
      {
        result.profile = envProfile;
        gzmsg << "[usv_hydro] Using profile from USV_HYDRO_PROFILE: "
              << result.profile << "\n";
      }
    }

    if (!configFile.empty())
    {
      std::string configError;
      if (!result.LoadFromFileProfile(configFile, result.profile, &configError))
      {
        gzerr << "[usv_hydro] Failed to load config profile: "
              << configError << "\n";
      }
    }

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
    if (_sdf->HasElement("scale_linear_drag_by_cell_count"))
    {
      result.scaleLinearDragByCellCount =
          _sdf->Get<bool>("scale_linear_drag_by_cell_count");
    }

    if (_sdf->HasElement("current_velocity"))
      result.currentVelocity = _sdf->Get<gz::math::Vector3d>("current_velocity");

    if (_sdf->HasElement("use_hydrostatic_stiffness_matrix"))
    {
      result.useHydrostaticStiffnessMatrix =
          _sdf->Get<bool>("use_hydrostatic_stiffness_matrix");
    }
    if (_sdf->HasElement("hydrostatic_stiffness_scale"))
    {
      result.hydrostaticStiffnessScale =
          _sdf->Get<double>("hydrostatic_stiffness_scale");
    }
    if (_sdf->HasElement("stiffness_heave_step"))
      result.stiffnessHeaveStep = _sdf->Get<double>("stiffness_heave_step");
    if (_sdf->HasElement("stiffness_angle_step"))
      result.stiffnessAngleStep = _sdf->Get<double>("stiffness_angle_step");

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

    gz::math::Vector3d effectiveLinearDrag = this->config.linearDrag;
    if (this->config.scaleLinearDragByCellCount && count > 0)
    {
      effectiveLinearDrag.X(effectiveLinearDrag.X() / count);
      effectiveLinearDrag.Y(effectiveLinearDrag.Y() / count);
      effectiveLinearDrag.Z(effectiveLinearDrag.Z() / count);
    }
    this->drag.SetParameters(
        this->config.fluidDensity, this->config.cd, effectiveLinearDrag);

    this->cellsReady = true;
    gzmsg << "[usv_hydro] Hydro cells initialized for link ["
          << this->config.linkName << "] count=" << count
          << " drag_scaled_by_cell_count="
          << (this->config.scaleLinearDragByCellCount ? "true" : "false")
          << " effective_linear_drag=[" << effectiveLinearDrag.X() << ", "
          << effectiveLinearDrag.Y() << ", " << effectiveLinearDrag.Z() << "]\n";
    return true;
  }

  private: bool configured{false};
  private: bool cellsReady{false};
  private: bool useHydrostaticStiffnessMatrix{false};
  private: bool hydrostaticStiffnessReady{false};

  private: double hydrostaticStiffnessScale{1.0};
  private: double stiffnessHeaveStep{0.01};
  private: double stiffnessAngleStep{1e-3};

  private: HydroConfig config;
  private: HydroCellGrid cellGrid;
  private: LinearHydrostaticStiffness hydrostaticStiffness;

  private: gz::math::Vector3d referencePositionWorld{0.0, 0.0, 0.0};
  private: gz::math::Quaterniond referenceRotationWorld;

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
