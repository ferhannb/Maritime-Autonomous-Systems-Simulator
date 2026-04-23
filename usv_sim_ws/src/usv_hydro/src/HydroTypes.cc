#include "usv_hydro/HydroTypes.hh"

#include <string>

namespace usv_hydro
{

bool HydroConfig::IsValid(std::string *_error) const
{
  if (this->linkName.empty())
  {
    if (_error)
      *_error = "link_name must not be empty";
    return false;
  }

  if (this->fluidDensity <= 0.0)
  {
    if (_error)
      *_error = "fluid_density must be > 0";
    return false;
  }

  if (this->gravity <= 0.0)
  {
    if (_error)
      *_error = "gravity must be > 0";
    return false;
  }

  if (this->cellsX <= 0 || this->cellsY <= 0 || this->cellsZ <= 0)
  {
    if (_error)
      *_error = "cells_x/cells_y/cells_z must be > 0";
    return false;
  }

  if (this->cd.X() < 0.0 || this->cd.Y() < 0.0 || this->cd.Z() < 0.0)
  {
    if (_error)
      *_error = "cd_x/cd_y/cd_z must be >= 0";
    return false;
  }

  if (this->linearDrag.X() < 0.0 ||
      this->linearDrag.Y() < 0.0 ||
      this->linearDrag.Z() < 0.0)
  {
    if (_error)
      *_error = "linear_drag_x/linear_drag_y/linear_drag_z must be >= 0";
    return false;
  }

  return true;
}

}  // namespace usv_hydro
