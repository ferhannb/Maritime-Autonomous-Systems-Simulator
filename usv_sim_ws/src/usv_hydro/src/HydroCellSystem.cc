#include <algorithm>
#include <array>
#include <chrono>
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
    // Row-major 6x6 hydrostatic restoring matrix.
    // DOF order: 0=surge, 1=sway, 2=heave, 3=roll, 4=pitch, 5=yaw
    std::array<double, 36> c{};
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
    this->useLinearSeakeepingModel = this->config.useLinearSeakeepingModel;
    this->seakeepingExcitationOmega = this->config.seakeepingExcitationOmega;
    this->seakeepingExcitationScale = this->config.seakeepingExcitationScale;
    this->seakeepingModelReady = false;

    if (this->useLinearSeakeepingModel)
    {
      std::string seakeepingError;
      if (!this->seakeepingModel.LoadFromFile(
              this->config.seakeepingCoeffsFile, &seakeepingError))
      {
        gzerr << "[usv_hydro] Failed to load linear seakeeping coeffs: "
              << seakeepingError << "\n";
        this->useLinearSeakeepingModel = false;
      }
      else
      {
        this->seakeepingModelReady = true;
        gzmsg << "[usv_hydro] Linear seakeeping model enabled "
              << "coeffs_file=[" << this->config.seakeepingCoeffsFile << "] "
              << "omega=" << this->seakeepingExcitationOmega
              << " excitation_scale=" << this->seakeepingExcitationScale
              << "\n";
      }
    }

    this->useCumminsRadiation = this->config.useCumminsRadiation;
    this->cumminsModelReady = false;
    if (this->useCumminsRadiation && !this->seakeepingModelReady)
    {
      gzerr << "[usv_hydro] Cummins radiation requires seakeeping coefficients "
            << "to be loaded successfully; disabling Cummins.\n";
      this->useCumminsRadiation = false;
    }
    else if (this->useCumminsRadiation)
    {
      gzmsg << "[usv_hydro] Cummins radiation model enabled "
            << "(will build on first PreUpdate step) "
            << "kernel_max_t=" << this->config.cumminsKernelMaxT
            << " kernel_dt=" << this->config.cumminsKernelDt << "\n";
    }

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
    this->previousBodyVelocity = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    this->previousBodyVelocityReady = false;
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
    const double dtSec = std::chrono::duration<double>(_info.dt).count();
    const double simTimeSec = std::chrono::duration<double>(_info.simTime).count();

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
              << "C33=" << this->hydrostaticStiffness.c[2 * 6 + 2]
              << " C44=" << this->hydrostaticStiffness.c[3 * 6 + 3]
              << " C55=" << this->hydrostaticStiffness.c[4 * 6 + 4] << "\n";
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

    // Build Cummins retardation kernel on first step (needs actual sim dt)
    if (this->useCumminsRadiation && !this->cumminsModelReady &&
        this->seakeepingModelReady && dtSec > 1e-9)
    {
      const double kdt = this->config.cumminsKernelDt > 0.0
          ? this->config.cumminsKernelDt
          : dtSec;
      std::string cumminsError;
      if (this->cumminsModel.BuildFromFrequencySamples(
              this->seakeepingModel.Samples(),
              this->config.cumminsKernelMaxT,
              kdt,
              &cumminsError))
      {
        this->cumminsModelReady = true;
        gzmsg << "[usv_hydro] Cummins radiation kernel built "
              << "kernel_dt=" << kdt
              << " kernel_max_t=" << this->config.cumminsKernelMaxT << "\n";
      }
      else
      {
        gzerr << "[usv_hydro] Failed to build Cummins radiation model: "
              << cumminsError << "; disabling.\n";
        this->useCumminsRadiation = false;
      }
    }

    // Body-frame kinematics — needed for seakeeping and/or Cummins
    const bool needsBodyKinematics =
        (this->useLinearSeakeepingModel && this->seakeepingModelReady) ||
        (this->useCumminsRadiation && this->cumminsModelReady);

    if (needsBodyKinematics)
    {
      const auto linearVelBody = pose->Rot().RotateVectorReverse(*linearVel);
      const auto angularVelBody = pose->Rot().RotateVectorReverse(*angularVel);

      const std::array<double, 6> bodyVelocity{
          linearVelBody.X(),
          linearVelBody.Y(),
          linearVelBody.Z(),
          angularVelBody.X(),
          angularVelBody.Y(),
          angularVelBody.Z()};

      std::array<double, 6> bodyAcceleration{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      if (this->previousBodyVelocityReady && dtSec > 1e-9)
      {
        for (std::size_t i = 0; i < bodyAcceleration.size(); ++i)
        {
          bodyAcceleration[i] =
              (bodyVelocity[i] - this->previousBodyVelocity[i]) / dtSec;
        }
      }
      this->previousBodyVelocity = bodyVelocity;
      this->previousBodyVelocityReady = true;

      if (this->useCumminsRadiation && this->cumminsModelReady)
      {
        // Cummins path: radiation via convolution + excitation separately
        const auto radWrench =
            this->cumminsModel.ComputeRadiationForce(bodyVelocity, bodyAcceleration);
        this->link.AddWorldWrench(
            _ecm,
            pose->Rot().RotateVector(
                gz::math::Vector3d(radWrench[0], radWrench[1], radWrench[2])),
            pose->Rot().RotateVector(
                gz::math::Vector3d(radWrench[3], radWrench[4], radWrench[5])));

        if (this->useLinearSeakeepingModel && this->seakeepingModelReady)
        {
          SeakeepingCoefficients coeffs;
          std::string evalError;
          if (this->seakeepingModel.Evaluate(
                  this->seakeepingExcitationOmega, &coeffs, &evalError))
          {
            const auto excWrench = this->seakeepingModel.ComputeExcitationWrench(
                coeffs, simTimeSec,
                this->seakeepingExcitationOmega,
                this->seakeepingExcitationScale);
            this->link.AddWorldWrench(
                _ecm,
                pose->Rot().RotateVector(
                    gz::math::Vector3d(excWrench[0], excWrench[1], excWrench[2])),
                pose->Rot().RotateVector(
                    gz::math::Vector3d(excWrench[3], excWrench[4], excWrench[5])));
          }
          else
          {
            gzerr << "[usv_hydro] Failed to evaluate seakeeping excitation: "
                  << evalError << "\n";
          }
        }
      }
      else if (this->useLinearSeakeepingModel && this->seakeepingModelReady)
      {
        // Legacy path: single-frequency radiation + excitation
        SeakeepingCoefficients coeffs;
        std::string evalError;
        if (this->seakeepingModel.Evaluate(
                this->seakeepingExcitationOmega, &coeffs, &evalError))
        {
          const auto bodyWrench = this->seakeepingModel.ComputeBodyWrench(
              coeffs,
              bodyVelocity,
              bodyAcceleration,
              simTimeSec,
              this->seakeepingExcitationOmega,
              this->seakeepingExcitationScale);
          this->link.AddWorldWrench(
              _ecm,
              pose->Rot().RotateVector(
                  gz::math::Vector3d(bodyWrench[0], bodyWrench[1], bodyWrench[2])),
              pose->Rot().RotateVector(
                  gz::math::Vector3d(bodyWrench[3], bodyWrench[4], bodyWrench[5])));
        }
        else
        {
          gzerr << "[usv_hydro] Failed to evaluate linear seakeeping coefficients: "
                << evalError << "\n";
        }
      }
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

    const auto aggregateAt = [&](
        const double _dx, const double _dy, const double _dz,
        const double _droll, const double _dpitch, const double _dyaw)
    {
      const auto positionWorld =
          this->referencePositionWorld +
          this->referenceRotationWorld.RotateVector(
              gz::math::Vector3d(_dx, _dy, _dz));
      const auto rotationWorld =
          this->referenceRotationWorld *
          gz::math::Quaterniond(_droll, _dpitch, _dyaw);
      const auto agg =
          this->ComputeHydrostaticAggregate(positionWorld, rotationWorld);
      // Express force and moment in reference body frame
      const auto fBody =
          this->referenceRotationWorld.RotateVectorReverse(agg.forceWorld);
      const auto mBody =
          this->referenceRotationWorld.RotateVectorReverse(agg.momentWorld);
      return std::array<double, 6>{
          fBody.X(), fBody.Y(), fBody.Z(),
          mBody.X(), mBody.Y(), mBody.Z()};
    };

    const double ds = this->stiffnessHeaveStep;   // translational perturbation
    const double da = this->stiffnessAngleStep;   // angular perturbation

    // For each column j, perturb DOF j positively and negatively.
    // DOF order: 0=surge, 1=sway, 2=heave, 3=roll, 4=pitch, 5=yaw
    struct ColPerturb
    {
      std::array<double, 6> plus;
      std::array<double, 6> minus;
      double step;
    };

    const std::array<ColPerturb, 6> perturb{{
      {aggregateAt( ds,  0,  0,  0,  0,  0), aggregateAt(-ds,  0,  0,  0,  0,  0), ds},
      {aggregateAt(  0, ds,  0,  0,  0,  0), aggregateAt(  0, -ds,  0,  0,  0,  0), ds},
      {aggregateAt(  0,  0, ds,  0,  0,  0), aggregateAt(  0,  0, -ds,  0,  0,  0), ds},
      {aggregateAt(  0,  0,  0, da,  0,  0), aggregateAt(  0,  0,  0, -da,  0,  0), da},
      {aggregateAt(  0,  0,  0,  0, da,  0), aggregateAt(  0,  0,  0,  0, -da,  0), da},
      {aggregateAt(  0,  0,  0,  0,  0, da), aggregateAt(  0,  0,  0,  0,  0, -da), da},
    }};

    auto &c = this->hydrostaticStiffness.c;
    c.fill(0.0);
    for (std::size_t col = 0; col < 6; ++col)
    {
      for (std::size_t row = 0; row < 6; ++row)
      {
        c[row * 6 + col] =
            -(perturb[col].plus[row] - perturb[col].minus[row]) /
            (2.0 * perturb[col].step);
      }
    }

    return std::isfinite(c[2 * 6 + 2]) &&   // C33: heave-heave
           std::isfinite(c[3 * 6 + 3]) &&   // C44: roll-roll
           std::isfinite(c[4 * 6 + 4]);     // C55: pitch-pitch
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

    // 6-DOF displacement vector in body (reference) frame
    // DOF order: 0=surge, 1=sway, 2=heave, 3=roll, 4=pitch, 5=yaw
    const std::array<double, 6> xi{
        deltaPositionRef.X(),
        deltaPositionRef.Y(),
        deltaPositionRef.Z(),
        eulerError.X(),
        eulerError.Y(),
        eulerError.Z()};

    const double scale = this->hydrostaticStiffnessScale;
    const auto &c = this->hydrostaticStiffness.c;

    std::array<double, 6> f{};
    f.fill(0.0);
    for (std::size_t row = 0; row < 6; ++row)
    {
      for (std::size_t col = 0; col < 6; ++col)
        f[row] -= scale * c[row * 6 + col] * xi[col];
    }

    // f[0..2] = force, f[3..5] = moment, in reference body frame
    const gz::math::Vector3d forceRef(f[0], f[1], f[2]);
    const gz::math::Vector3d torqueRef(f[3], f[4], f[5]);

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

    if (_sdf->HasElement("use_linear_seakeeping_model"))
    {
      result.useLinearSeakeepingModel =
          _sdf->Get<bool>("use_linear_seakeeping_model");
    }
    if (_sdf->HasElement("seakeeping_coeffs_file"))
    {
      result.seakeepingCoeffsFile =
          _sdf->Get<std::string>("seakeeping_coeffs_file");
    }
    if (_sdf->HasElement("seakeeping_excitation_omega"))
    {
      result.seakeepingExcitationOmega =
          _sdf->Get<double>("seakeeping_excitation_omega");
    }
    if (_sdf->HasElement("seakeeping_excitation_scale"))
    {
      result.seakeepingExcitationScale =
          _sdf->Get<double>("seakeeping_excitation_scale");
    }
    if (_sdf->HasElement("cummins_radiation_enabled"))
      result.useCumminsRadiation = _sdf->Get<bool>("cummins_radiation_enabled");
    if (_sdf->HasElement("cummins_kernel_max_t"))
      result.cumminsKernelMaxT = _sdf->Get<double>("cummins_kernel_max_t");
    if (_sdf->HasElement("cummins_kernel_dt"))
      result.cumminsKernelDt = _sdf->Get<double>("cummins_kernel_dt");

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
  private: bool useLinearSeakeepingModel{false};
  private: bool seakeepingModelReady{false};
  private: bool useCumminsRadiation{false};
  private: bool cumminsModelReady{false};
  private: bool previousBodyVelocityReady{false};

  private: double hydrostaticStiffnessScale{1.0};
  private: double stiffnessHeaveStep{0.01};
  private: double stiffnessAngleStep{1e-3};
  private: double seakeepingExcitationOmega{0.0};
  private: double seakeepingExcitationScale{1.0};

  private: HydroConfig config;
  private: HydroCellGrid cellGrid;
  private: LinearHydrostaticStiffness hydrostaticStiffness;
  private: LinearSeakeepingModel seakeepingModel;
  private: CumminsRadiationModel cumminsModel;
  private: std::array<double, 6> previousBodyVelocity{
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

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
