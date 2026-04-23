#include "usv_hydro/DragModel.hh"

#include <cmath>

namespace usv_hydro
{

DragModel::DragModel(
    const double _fluidDensity,
    const gz::math::Vector3d &_cd,
    const gz::math::Vector3d &_linearDrag)
{
  this->SetParameters(_fluidDensity, _cd, _linearDrag);
}

void DragModel::SetParameters(
    const double _fluidDensity,
    const gz::math::Vector3d &_cd,
    const gz::math::Vector3d &_linearDrag)
{
  this->fluidDensity = _fluidDensity;
  this->cd = _cd;
  this->linearDrag = _linearDrag;
}

gz::math::Vector3d DragModel::ComputeForceBody(
    const gz::math::Vector3d &_relativeVelocityBody,
    const gz::math::Vector3d &_dragAreaCell) const
{
  gz::math::Vector3d dragBody;

  dragBody.X(
      -this->linearDrag.X() * _relativeVelocityBody.X()
      - 0.5 * this->fluidDensity * this->cd.X() * _dragAreaCell.X()
      * std::abs(_relativeVelocityBody.X()) * _relativeVelocityBody.X());

  dragBody.Y(
      -this->linearDrag.Y() * _relativeVelocityBody.Y()
      - 0.5 * this->fluidDensity * this->cd.Y() * _dragAreaCell.Y()
      * std::abs(_relativeVelocityBody.Y()) * _relativeVelocityBody.Y());

  dragBody.Z(
      -this->linearDrag.Z() * _relativeVelocityBody.Z()
      - 0.5 * this->fluidDensity * this->cd.Z() * _dragAreaCell.Z()
      * std::abs(_relativeVelocityBody.Z()) * _relativeVelocityBody.Z());

  return dragBody;
}

}  // namespace usv_hydro
