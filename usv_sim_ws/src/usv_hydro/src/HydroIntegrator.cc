#include "usv_hydro/HydroIntegrator.hh"

namespace usv_hydro
{

HydroCellForces HydroIntegrator::ComputeCellForces(
    const gz::math::Vector3d &_offsetLocal,
    const HydroCellGrid &_grid,
    const HydroKinematics &_kinematics,
    const EnvironmentModel &_environment,
    const BuoyancyModel &_buoyancy,
    const DragModel &_drag) const
{
  HydroCellForces result;

  const auto offsetWorld = _kinematics.rotationWorld.RotateVector(_offsetLocal);
  const auto cellPositionWorld = _kinematics.positionWorld + offsetWorld;

  const double depth = _environment.DepthAt(cellPositionWorld.Z());
  const double submergence =
      _buoyancy.ComputeSubmergence(depth, _grid.cellHeightApprox);
  if (submergence <= 0.0)
    return result;

  result.submerged = true;
  result.submergence = submergence;
  result.buoyancyWorld =
      _buoyancy.ComputeForceWorld(_grid.cellVolume, submergence);

  const auto cellVelocityWorld =
      _kinematics.linearVelocityWorld
      + _kinematics.angularVelocityWorld.Cross(offsetWorld);
  const auto relativeVelocityWorld =
      _environment.RelativeVelocityWorld(cellVelocityWorld);
  const auto relativeVelocityBody =
      _kinematics.rotationWorld.RotateVectorReverse(relativeVelocityWorld);

  const auto dragBody =
      _drag.ComputeForceBody(relativeVelocityBody, _grid.dragAreaCell);
  result.dragWorld = _kinematics.rotationWorld.RotateVector(
      dragBody * submergence);

  return result;
}

}  // namespace usv_hydro
