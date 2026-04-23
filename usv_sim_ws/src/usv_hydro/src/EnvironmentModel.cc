#include "usv_hydro/EnvironmentModel.hh"

namespace usv_hydro
{

EnvironmentModel::EnvironmentModel(const HydroConfig &_config)
{
  this->SetFromConfig(_config);
}

void EnvironmentModel::SetFromConfig(const HydroConfig &_config)
{
  this->waterLevel = _config.waterLevel;
  this->fluidDensity = _config.fluidDensity;
  this->gravity = _config.gravity;
  this->currentVelocity = _config.currentVelocity;
}

double EnvironmentModel::DepthAt(const double _worldZ) const
{
  return this->waterLevel - _worldZ;
}

gz::math::Vector3d EnvironmentModel::RelativeVelocityWorld(
    const gz::math::Vector3d &_velocityWorld) const
{
  return _velocityWorld - this->currentVelocity;
}

double EnvironmentModel::FluidDensity() const
{
  return this->fluidDensity;
}

double EnvironmentModel::Gravity() const
{
  return this->gravity;
}

}  // namespace usv_hydro
