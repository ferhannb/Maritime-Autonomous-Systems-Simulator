#ifndef USV_HYDRO_BUOYANCY_MODEL_HH_
#define USV_HYDRO_BUOYANCY_MODEL_HH_

#include <gz/math/Vector3.hh>

namespace usv_hydro
{

class BuoyancyModel
{
  public: BuoyancyModel() = default;
  public: BuoyancyModel(const double _fluidDensity, const double _gravity);

  public: void SetParameters(const double _fluidDensity, const double _gravity);

  public: double ComputeSubmergence(
      const double _depth, const double _cellHeight) const;
  public: gz::math::Vector3d ComputeForceWorld(
      const double _cellVolume, const double _submergence) const;

  private: double fluidDensity{1025.0};
  private: double gravity{9.81};
};

}  // namespace usv_hydro

#endif  // USV_HYDRO_BUOYANCY_MODEL_HH_
