#ifndef USV_HYDRO_HYDRO_INTEGRATOR_HH_
#define USV_HYDRO_HYDRO_INTEGRATOR_HH_

#include <gz/math/Quaternion.hh>
#include <gz/math/Vector3.hh>

#include "usv_hydro/BuoyancyModel.hh"
#include "usv_hydro/DragModel.hh"
#include "usv_hydro/EnvironmentModel.hh"
#include "usv_hydro/HydroTypes.hh"

namespace usv_hydro
{

struct HydroKinematics
{
  gz::math::Vector3d positionWorld{0.0, 0.0, 0.0};
  gz::math::Quaterniond rotationWorld;
  gz::math::Vector3d linearVelocityWorld{0.0, 0.0, 0.0};
  gz::math::Vector3d angularVelocityWorld{0.0, 0.0, 0.0};
};

class HydroIntegrator
{
  public: HydroCellForces ComputeCellForces(
      const gz::math::Vector3d &_offsetLocal,
      const HydroCellGrid &_grid,
      const HydroKinematics &_kinematics,
      const EnvironmentModel &_environment,
      const BuoyancyModel &_buoyancy,
      const DragModel &_drag) const;
};

}  // namespace usv_hydro

#endif  // USV_HYDRO_HYDRO_INTEGRATOR_HH_
