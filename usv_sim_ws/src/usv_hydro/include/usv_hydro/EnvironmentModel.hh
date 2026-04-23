#ifndef USV_HYDRO_ENVIRONMENT_MODEL_HH_
#define USV_HYDRO_ENVIRONMENT_MODEL_HH_

#include <gz/math/Vector3.hh>

#include "usv_hydro/HydroTypes.hh"

namespace usv_hydro
{

class EnvironmentModel
{
  public: EnvironmentModel() = default;
  public: explicit EnvironmentModel(const HydroConfig &_config);

  public: void SetFromConfig(const HydroConfig &_config);

  public: double DepthAt(const double _worldZ) const;
  public: gz::math::Vector3d RelativeVelocityWorld(
      const gz::math::Vector3d &_velocityWorld) const;

  public: double FluidDensity() const;
  public: double Gravity() const;

  private: double waterLevel{0.0};
  private: double fluidDensity{1000.0};
  private: double gravity{9.81};
  private: gz::math::Vector3d currentVelocity{0.0, 0.0, 0.0};
};

}  // namespace usv_hydro

#endif  // USV_HYDRO_ENVIRONMENT_MODEL_HH_
