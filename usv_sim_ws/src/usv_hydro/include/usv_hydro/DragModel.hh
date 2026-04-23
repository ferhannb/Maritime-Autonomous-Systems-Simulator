#ifndef USV_HYDRO_DRAG_MODEL_HH_
#define USV_HYDRO_DRAG_MODEL_HH_

#include <gz/math/Vector3.hh>

namespace usv_hydro
{

class DragModel
{
  public: DragModel() = default;
  public: DragModel(
      const double _fluidDensity,
      const gz::math::Vector3d &_cd,
      const gz::math::Vector3d &_linearDrag);

  public: void SetParameters(
      const double _fluidDensity,
      const gz::math::Vector3d &_cd,
      const gz::math::Vector3d &_linearDrag);

  public: gz::math::Vector3d ComputeForceBody(
      const gz::math::Vector3d &_relativeVelocityBody,
      const gz::math::Vector3d &_dragAreaCell) const;

  private: double fluidDensity{1000.0};
  private: gz::math::Vector3d cd{0.9, 1.2, 2.0};
  private: gz::math::Vector3d linearDrag{20.0, 40.0, 70.0};
};

}  // namespace usv_hydro

#endif  // USV_HYDRO_DRAG_MODEL_HH_
