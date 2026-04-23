#include "usv_hydro/BuoyancyModel.hh"

#include <algorithm>

namespace usv_hydro
{

BuoyancyModel::BuoyancyModel(const double _fluidDensity, const double _gravity)
{
  this->SetParameters(_fluidDensity, _gravity);
}

void BuoyancyModel::SetParameters(const double _fluidDensity, const double _gravity)
{
  this->fluidDensity = _fluidDensity;
  this->gravity = _gravity;
}

double BuoyancyModel::ComputeSubmergence(
    const double _depth, const double _cellHeight) const
{
  if (_cellHeight <= 1e-9)
    return 0.0;

  // `_depth` is measured at the cell center. Convert to a fractional fill
  // of the full cell height: 0 at top face touching SWL, 1 at bottom face.
  return std::clamp(_depth / _cellHeight + 0.5, 0.0, 1.0);
}

gz::math::Vector3d BuoyancyModel::ComputeForceWorld(
    const double _cellVolume, const double _submergence) const
{
  const double displacedVolume = _cellVolume * _submergence;
  return gz::math::Vector3d(
      0.0, 0.0, this->fluidDensity * this->gravity * displacedVolume);
}

}  // namespace usv_hydro
