#include "usv_hydro/HydroTypes.hh"

#include <fstream>
#include <string>

#include <yaml-cpp/yaml.h>

namespace usv_hydro
{

namespace
{

bool ParseVector3(
    const YAML::Node &_node,
    gz::math::Vector3d *_out,
    const std::string &_field,
    std::string *_error)
{
  if (!_node || !_node.IsSequence() || _node.size() != 3)
  {
    if (_error)
      *_error = _field + " must be a sequence of 3 numbers";
    return false;
  }

  try
  {
    _out->Set(
        _node[0].as<double>(),
        _node[1].as<double>(),
        _node[2].as<double>());
  }
  catch (const YAML::Exception &)
  {
    if (_error)
      *_error = _field + " contains invalid numeric values";
    return false;
  }
  return true;
}

bool ParseCells(
    const YAML::Node &_node,
    int *_x,
    int *_y,
    int *_z,
    std::string *_error)
{
  if (!_node || !_node.IsSequence() || _node.size() != 3)
  {
    if (_error)
      *_error = "cells must be a sequence of 3 integers";
    return false;
  }

  try
  {
    *_x = _node[0].as<int>();
    *_y = _node[1].as<int>();
    *_z = _node[2].as<int>();
  }
  catch (const YAML::Exception &)
  {
    if (_error)
      *_error = "cells contains invalid integer values";
    return false;
  }
  return true;
}

bool OverlayConfigNode(
    const YAML::Node &_node,
    HydroConfig *_cfg,
    std::string *_error)
{
  if (!_node)
    return true;
  if (!_node.IsMap())
  {
    if (_error)
      *_error = "config section must be a map";
    return false;
  }

  try
  {
    if (_node["fluid_density"])
      _cfg->fluidDensity = _node["fluid_density"].as<double>();
    if (_node["water_level"])
      _cfg->waterLevel = _node["water_level"].as<double>();
    if (_node["gravity"])
      _cfg->gravity = _node["gravity"].as<double>();
    if (_node["link_name"])
      _cfg->linkName = _node["link_name"].as<std::string>();

    if (_node["cells"] &&
        !ParseCells(_node["cells"], &_cfg->cellsX, &_cfg->cellsY, &_cfg->cellsZ, _error))
      return false;

    if (_node["cd"] && !ParseVector3(_node["cd"], &_cfg->cd, "cd", _error))
      return false;

    if (_node["current_velocity"] &&
        !ParseVector3(
            _node["current_velocity"], &_cfg->currentVelocity, "current_velocity", _error))
      return false;

    if (_node["drag"])
    {
      const auto drag = _node["drag"];
      if (!drag.IsMap())
      {
        if (_error)
          *_error = "drag must be a map";
        return false;
      }
      if (drag["linear_total"] &&
          !ParseVector3(drag["linear_total"], &_cfg->linearDrag, "drag.linear_total", _error))
        return false;
      if (drag["scale_by_cell_count"])
        _cfg->scaleLinearDragByCellCount = drag["scale_by_cell_count"].as<bool>();
    }

    if (_node["hydrostatic_stiffness"])
    {
      const auto stiffness = _node["hydrostatic_stiffness"];
      if (!stiffness.IsMap())
      {
        if (_error)
          *_error = "hydrostatic_stiffness must be a map";
        return false;
      }
      if (stiffness["enabled"])
        _cfg->useHydrostaticStiffnessMatrix = stiffness["enabled"].as<bool>();
      if (stiffness["scale"])
        _cfg->hydrostaticStiffnessScale = stiffness["scale"].as<double>();
      if (stiffness["heave_step"])
        _cfg->stiffnessHeaveStep = stiffness["heave_step"].as<double>();
      if (stiffness["angle_step"])
        _cfg->stiffnessAngleStep = stiffness["angle_step"].as<double>();
    }
  }
  catch (const YAML::Exception &ex)
  {
    if (_error)
      *_error = std::string("yaml parse error: ") + ex.what();
    return false;
  }

  return true;
}

}  // namespace

bool HydroConfig::LoadFromFileProfile(
    const std::string &_configFile,
    const std::string &_profile,
    std::string *_error)
{
  if (_configFile.empty())
  {
    if (_error)
      *_error = "config_file must not be empty";
    return false;
  }

  std::ifstream in(_configFile);
  if (!in.good())
  {
    if (_error)
      *_error = "could not open config_file: " + _configFile;
    return false;
  }

  YAML::Node root;
  try
  {
    root = YAML::Load(in);
  }
  catch (const YAML::Exception &ex)
  {
    if (_error)
      *_error = std::string("failed to load yaml: ") + ex.what();
    return false;
  }

  if (!root || !root.IsMap())
  {
    if (_error)
      *_error = "yaml root must be a map";
    return false;
  }

  if (!OverlayConfigNode(root["defaults"], this, _error))
    return false;

  const std::string selectedProfile = _profile.empty() ? "baseline" : _profile;
  if (root["profiles"])
  {
    if (!root["profiles"].IsMap())
    {
      if (_error)
        *_error = "profiles must be a map";
      return false;
    }

    if (root["profiles"][selectedProfile])
    {
      if (!OverlayConfigNode(root["profiles"][selectedProfile], this, _error))
        return false;
    }
    else
    {
      if (_error)
        *_error = "profile not found in config file: " + selectedProfile;
      return false;
    }
  }
  else if (selectedProfile != "baseline")
  {
    if (_error)
      *_error = "profiles section missing in config file";
    return false;
  }

  this->configFile = _configFile;
  this->profile = selectedProfile;
  return true;
}

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

  if (this->hydrostaticStiffnessScale < 0.0)
  {
    if (_error)
      *_error = "hydrostatic_stiffness_scale must be >= 0";
    return false;
  }

  if (this->stiffnessHeaveStep <= 0.0)
  {
    if (_error)
      *_error = "stiffness_heave_step must be > 0";
    return false;
  }

  if (this->stiffnessAngleStep <= 0.0)
  {
    if (_error)
      *_error = "stiffness_angle_step must be > 0";
    return false;
  }

  return true;
}

}  // namespace usv_hydro
