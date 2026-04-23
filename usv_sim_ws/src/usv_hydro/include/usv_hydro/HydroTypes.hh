#ifndef USV_HYDRO_HYDRO_TYPES_HH_
#define USV_HYDRO_HYDRO_TYPES_HH_

#include <string>
#include <vector>
#include <gz/math/Vector3.hh>

namespace usv_hydro
{

struct HydroConfig
{
  std::string linkName{"hull_link"};

  std::string configFile{};
  std::string profile{"baseline"};

  double fluidDensity{1000.0};
  double waterLevel{0.0};
  double gravity{9.81};

  int cellsX{12};
  int cellsY{6};
  int cellsZ{4};

  gz::math::Vector3d cd{0.9, 1.2, 2.0};
  gz::math::Vector3d linearDrag{20.0, 40.0, 70.0};
  bool scaleLinearDragByCellCount{true};
  gz::math::Vector3d currentVelocity{0.0, 0.0, 0.0};

  bool useHydrostaticStiffnessMatrix{false};
  double hydrostaticStiffnessScale{1.0};
  double stiffnessHeaveStep{0.01};
  double stiffnessAngleStep{1e-3};

  bool LoadFromFileProfile(
      const std::string &_configFile,
      const std::string &_profile,
      std::string *_error = nullptr);

  bool IsValid(std::string *_error = nullptr) const;
};

struct HydroCellGrid
{
  std::vector<gz::math::Vector3d> offsets;
  double cellVolume{0.0};
  double cellHeightApprox{0.1};
  gz::math::Vector3d dragAreaCell{0.0, 0.0, 0.0};
};

struct HydroCellForces
{
  gz::math::Vector3d buoyancyWorld{0.0, 0.0, 0.0};
  gz::math::Vector3d dragWorld{0.0, 0.0, 0.0};
  double submergence{0.0};
  bool submerged{false};
};

}  // namespace usv_hydro

#endif  // USV_HYDRO_HYDRO_TYPES_HH_
