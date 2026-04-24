#include "usv_hydro/HydroTypes.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <string>
#include <utility>

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

bool ParseVector6(
    const YAML::Node &_node,
    std::array<double, 6> *_out,
    const std::string &_field,
    std::string *_error)
{
  if (!_node || !_node.IsSequence() || _node.size() != 6)
  {
    if (_error)
      *_error = _field + " must be a sequence of 6 numbers";
    return false;
  }

  try
  {
    for (std::size_t i = 0; i < 6; ++i)
      (*_out)[i] = _node[i].as<double>();
  }
  catch (const YAML::Exception &)
  {
    if (_error)
      *_error = _field + " contains invalid numeric values";
    return false;
  }
  return true;
}

bool ParseMatrix6x6(
    const YAML::Node &_node,
    std::array<double, 36> *_out,
    const std::string &_field,
    std::string *_error)
{
  if (!_node)
    return false;

  try
  {
    if (_node.IsSequence() && _node.size() == 36)
    {
      for (std::size_t i = 0; i < 36; ++i)
        (*_out)[i] = _node[i].as<double>();
      return true;
    }

    if (_node.IsSequence() && _node.size() == 6)
    {
      for (std::size_t r = 0; r < 6; ++r)
      {
        const auto row = _node[r];
        if (!row.IsSequence() || row.size() != 6)
        {
          if (_error)
            *_error = _field + " row must contain 6 numbers";
          return false;
        }
        for (std::size_t c = 0; c < 6; ++c)
          (*_out)[r * 6 + c] = row[c].as<double>();
      }
      return true;
    }
  }
  catch (const YAML::Exception &)
  {
    if (_error)
      *_error = _field + " contains invalid numeric values";
    return false;
  }

  if (_error)
    *_error = _field + " must be a 36-value flat sequence or 6x6 nested sequence";
  return false;
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

    if (_node["linear_seakeeping"])
    {
      const auto seakeeping = _node["linear_seakeeping"];
      if (!seakeeping.IsMap())
      {
        if (_error)
          *_error = "linear_seakeeping must be a map";
        return false;
      }
      if (seakeeping["enabled"])
        _cfg->useLinearSeakeepingModel = seakeeping["enabled"].as<bool>();
      if (seakeeping["coeffs_file"])
        _cfg->seakeepingCoeffsFile = seakeeping["coeffs_file"].as<std::string>();
      if (seakeeping["excitation_omega"])
      {
        _cfg->seakeepingExcitationOmega =
            seakeeping["excitation_omega"].as<double>();
      }
      if (seakeeping["excitation_scale"])
      {
        _cfg->seakeepingExcitationScale =
            seakeeping["excitation_scale"].as<double>();
      }
    }

    if (_node["cummins_radiation"])
    {
      const auto cummins = _node["cummins_radiation"];
      if (!cummins.IsMap())
      {
        if (_error)
          *_error = "cummins_radiation must be a map";
        return false;
      }
      if (cummins["enabled"])
        _cfg->useCumminsRadiation = cummins["enabled"].as<bool>();
      if (cummins["kernel_max_t"])
        _cfg->cumminsKernelMaxT = cummins["kernel_max_t"].as<double>();
      if (cummins["kernel_dt"])
        _cfg->cumminsKernelDt = cummins["kernel_dt"].as<double>();
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

  if (!this->seakeepingCoeffsFile.empty())
  {
    const std::filesystem::path coeffPath(this->seakeepingCoeffsFile);
    if (coeffPath.is_relative())
    {
      const auto baseDir = std::filesystem::path(_configFile).parent_path();
      this->seakeepingCoeffsFile = (baseDir / coeffPath).lexically_normal().string();
    }
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

  if (this->useLinearSeakeepingModel && this->seakeepingCoeffsFile.empty())
  {
    if (_error)
      *_error = "seakeeping_coeffs_file must not be empty when linear_seakeeping is enabled";
    return false;
  }

  if (this->seakeepingExcitationOmega < 0.0)
  {
    if (_error)
      *_error = "seakeeping_excitation_omega must be >= 0";
    return false;
  }

  if (this->seakeepingExcitationScale < 0.0)
  {
    if (_error)
      *_error = "seakeeping_excitation_scale must be >= 0";
    return false;
  }

  if (this->useCumminsRadiation)
  {
    if (!this->useLinearSeakeepingModel)
    {
      if (_error)
        *_error =
            "cummins_radiation requires linear_seakeeping to be enabled "
            "(needs the B(omega) frequency samples)";
      return false;
    }

    if (this->cumminsKernelMaxT <= 0.0)
    {
      if (_error)
        *_error = "cummins_kernel_max_t must be > 0";
      return false;
    }

    if (this->cumminsKernelDt < 0.0)
    {
      if (_error)
        *_error = "cummins_kernel_dt must be >= 0 (0 = use simulation dt)";
      return false;
    }
  }

  return true;
}

bool LinearSeakeepingModel::LoadFromFile(
    const std::string &_coeffFile,
    std::string *_error)
{
  if (_coeffFile.empty())
  {
    if (_error)
      *_error = "seakeeping coeff file path is empty";
    return false;
  }

  std::ifstream in(_coeffFile);
  if (!in.good())
  {
    if (_error)
      *_error = "could not open seakeeping coeff file: " + _coeffFile;
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
      *_error = std::string("failed to load seakeeping yaml: ") + ex.what();
    return false;
  }

  const auto frequencies = root["frequencies"];
  if (!frequencies || !frequencies.IsSequence() || frequencies.size() == 0)
  {
    if (_error)
      *_error = "frequencies must be a non-empty sequence";
    return false;
  }

  std::vector<FrequencySample> parsed;
  parsed.reserve(frequencies.size());

  for (std::size_t i = 0; i < frequencies.size(); ++i)
  {
    const auto node = frequencies[i];
    if (!node.IsMap())
    {
      if (_error)
        *_error = "each frequencies entry must be a map";
      return false;
    }
    if (!node["omega"])
    {
      if (_error)
        *_error = "each frequencies entry must define omega";
      return false;
    }

    FrequencySample sample;
    try
    {
      sample.omega = node["omega"].as<double>();
    }
    catch (const YAML::Exception &)
    {
      if (_error)
        *_error = "omega contains invalid numeric value";
      return false;
    }
    if (sample.omega < 0.0)
    {
      if (_error)
        *_error = "omega must be >= 0";
      return false;
    }

    if (!ParseMatrix6x6(
            node["added_mass"],
            &sample.coeffs.addedMass,
            "frequencies.added_mass",
            _error))
      return false;
    if (!ParseMatrix6x6(
            node["damping"],
            &sample.coeffs.damping,
            "frequencies.damping",
            _error))
      return false;

    if (node["excitation_re"] &&
        !ParseVector6(
            node["excitation_re"],
            &sample.coeffs.excitationRe,
            "frequencies.excitation_re",
            _error))
      return false;
    if (node["excitation_im"] &&
        !ParseVector6(
            node["excitation_im"],
            &sample.coeffs.excitationIm,
            "frequencies.excitation_im",
            _error))
      return false;

    parsed.emplace_back(std::move(sample));
  }

  std::sort(
      parsed.begin(),
      parsed.end(),
      [](const FrequencySample &_a, const FrequencySample &_b)
      {
        return _a.omega < _b.omega;
      });

  for (std::size_t i = 1; i < parsed.size(); ++i)
  {
    if (std::abs(parsed[i].omega - parsed[i - 1].omega) < 1e-12)
    {
      if (_error)
        *_error = "frequencies contain duplicate omega entries";
      return false;
    }
  }

  this->samples = std::move(parsed);
  return true;
}

bool LinearSeakeepingModel::IsLoaded() const
{
  return !this->samples.empty();
}

bool LinearSeakeepingModel::Evaluate(
    const double _omega,
    SeakeepingCoefficients *_coeffs,
    std::string *_error) const
{
  if (!_coeffs)
  {
    if (_error)
      *_error = "output coeffs pointer is null";
    return false;
  }

  if (this->samples.empty())
  {
    if (_error)
      *_error = "seakeeping model has no loaded samples";
    return false;
  }

  if (this->samples.size() == 1 || _omega <= this->samples.front().omega)
  {
    *_coeffs = this->samples.front().coeffs;
    return true;
  }

  if (_omega >= this->samples.back().omega)
  {
    *_coeffs = this->samples.back().coeffs;
    return true;
  }

  const auto upperIt = std::upper_bound(
      this->samples.begin(),
      this->samples.end(),
      _omega,
      [](const double _value, const FrequencySample &_sample)
      {
        return _value < _sample.omega;
      });

  const std::size_t hi = static_cast<std::size_t>(
      std::distance(this->samples.begin(), upperIt));
  const std::size_t lo = hi - 1;

  const double omegaLo = this->samples[lo].omega;
  const double omegaHi = this->samples[hi].omega;
  const double denom = omegaHi - omegaLo;
  if (denom <= 1e-12)
  {
    *_coeffs = this->samples[lo].coeffs;
    return true;
  }

  const double alpha = (_omega - omegaLo) / denom;
  const auto &cLo = this->samples[lo].coeffs;
  const auto &cHi = this->samples[hi].coeffs;

  for (std::size_t i = 0; i < 36; ++i)
  {
    _coeffs->addedMass[i] =
        (1.0 - alpha) * cLo.addedMass[i] + alpha * cHi.addedMass[i];
    _coeffs->damping[i] =
        (1.0 - alpha) * cLo.damping[i] + alpha * cHi.damping[i];
  }

  for (std::size_t i = 0; i < 6; ++i)
  {
    _coeffs->excitationRe[i] =
        (1.0 - alpha) * cLo.excitationRe[i] + alpha * cHi.excitationRe[i];
    _coeffs->excitationIm[i] =
        (1.0 - alpha) * cLo.excitationIm[i] + alpha * cHi.excitationIm[i];
  }

  return true;
}

  std::array<double, 6> LinearSeakeepingModel::ComputeBodyWrench(
    const SeakeepingCoefficients &_coeffs,
    const std::array<double, 6> &_velocityBody,
    const std::array<double, 6> &_accelerationBody,
    const double _timeSec,
    const double _omega,
    const double _excitationScale) const
{
  std::array<double, 6> wrench{};
  const double c = std::cos(_omega * _timeSec);
  const double s = std::sin(_omega * _timeSec);

  for (std::size_t row = 0; row < 6; ++row)
  {
    const std::size_t rowBase = row * 6;
    double hydroTerm = 0.0;
    for (std::size_t col = 0; col < 6; ++col)
    {
      hydroTerm +=
          _coeffs.addedMass[rowBase + col] * _accelerationBody[col];
      hydroTerm += _coeffs.damping[rowBase + col] * _velocityBody[col];
    }

    const double excitation = _excitationScale *
        (_coeffs.excitationRe[row] * c - _coeffs.excitationIm[row] * s);
    wrench[row] = -hydroTerm + excitation;
  }

  return wrench;
}

std::array<double, 6> LinearSeakeepingModel::ComputeExcitationWrench(
    const SeakeepingCoefficients &_coeffs,
    const double _timeSec,
    const double _omega,
    const double _excitationScale) const
{
  std::array<double, 6> wrench{};
  const double c = std::cos(_omega * _timeSec);
  const double s = std::sin(_omega * _timeSec);

  for (std::size_t row = 0; row < 6; ++row)
  {
    wrench[row] = _excitationScale *
        (_coeffs.excitationRe[row] * c - _coeffs.excitationIm[row] * s);
  }

  return wrench;
}

const std::vector<LinearSeakeepingModel::FrequencySample> &
LinearSeakeepingModel::Samples() const
{
  return this->samples;
}

// ---------------------------------------------------------------------------
// CumminsRadiationModel
// ---------------------------------------------------------------------------

bool CumminsRadiationModel::BuildFromFrequencySamples(
    const std::vector<LinearSeakeepingModel::FrequencySample> &_samples,
    const double _maxT,
    const double _kernelDt,
    std::string *_error)
{
  if (_samples.size() < 2)
  {
    if (_error)
      *_error =
          "at least 2 frequency samples required to compute retardation kernel";
    return false;
  }

  if (_maxT <= 0.0)
  {
    if (_error)
      *_error = "kernel_max_t must be > 0";
    return false;
  }

  if (_kernelDt <= 0.0)
  {
    if (_error)
      *_error = "kernel_dt must be > 0";
    return false;
  }

  this->kernelDt = _kernelDt;
  this->maxHistoryLen =
      static_cast<std::size_t>(std::ceil(_maxT / _kernelDt)) + 1u;

  const std::size_t nT = this->maxHistoryLen;
  const std::size_t nOmega = _samples.size();

  // K_ij(t) = (2/pi) * trapz[ B_ij(omega) * cos(omega * t), omega ]
  // Using proper non-uniform trapezoidal rule over the omega axis.
  this->kernel.assign(nT, {});
  for (std::size_t k = 0; k < nT; ++k)
  {
    const double t = static_cast<double>(k) * _kernelDt;
    auto &Kt = this->kernel[k];
    Kt.fill(0.0);

    for (std::size_t n = 0; n + 1 < nOmega; ++n)
    {
      const double omLo = _samples[n].omega;
      const double omHi = _samples[n + 1].omega;
      const double dOm  = omHi - omLo;
      const double cosLo = std::cos(omLo * t);
      const double cosHi = std::cos(omHi * t);

      for (std::size_t ij = 0; ij < 36; ++ij)
      {
        Kt[ij] += 0.5 * dOm *
            (_samples[n].coeffs.damping[ij] * cosLo +
             _samples[n + 1].coeffs.damping[ij] * cosHi);
      }
    }

    for (std::size_t ij = 0; ij < 36; ++ij)
      Kt[ij] *= (2.0 / M_PI);
  }

  // Truncate kernel: drop tail where all elements are smaller than 1e-6 * K(0)
  if (nT > 3)
  {
    double maxK0 = 0.0;
    for (std::size_t ij = 0; ij < 36; ++ij)
      maxK0 = std::max(maxK0, std::abs(this->kernel[0][ij]));

    if (maxK0 > 0.0)
    {
      std::size_t lastSig = nT - 1;
      for (std::size_t k = nT - 1; k > 2; --k)
      {
        double maxKt = 0.0;
        for (std::size_t ij = 0; ij < 36; ++ij)
          maxKt = std::max(maxKt, std::abs(this->kernel[k][ij]));
        if (maxKt > 1e-6 * maxK0)
        {
          lastSig = k;
          break;
        }
      }
      if (lastSig < nT - 1)
      {
        this->kernel.resize(lastSig + 1);
        this->maxHistoryLen = lastSig + 1;
      }
    }
  }

  // Compute A(inf) via the Ogilvie relation:
  //   A_inf = A(omega_max) + (1/omega_max) * integral_0^T K(t)*sin(omega_max*t) dt
  const double omegaMax = _samples.back().omega;
  this->addedMassInf = _samples.back().coeffs.addedMass;

  if (omegaMax > 1e-12 && !this->kernel.empty())
  {
    const std::size_t nK = this->kernel.size();
    std::array<double, 36> integralKSin{};
    integralKSin.fill(0.0);

    for (std::size_t k = 0; k < nK; ++k)
    {
      const double t = static_cast<double>(k) * _kernelDt;
      const double sinTerm = std::sin(omegaMax * t);
      const double w = (k == 0 || k == nK - 1) ? 0.5 : 1.0;
      for (std::size_t ij = 0; ij < 36; ++ij)
        integralKSin[ij] += w * this->kernel[k][ij] * sinTerm;
    }

    for (std::size_t ij = 0; ij < 36; ++ij)
    {
      integralKSin[ij] *= _kernelDt;
      this->addedMassInf[ij] += integralKSin[ij] / omegaMax;
    }
  }

  this->velHistory.clear();
  this->ready = true;
  return true;
}

bool CumminsRadiationModel::IsReady() const
{
  return this->ready;
}

std::array<double, 6> CumminsRadiationModel::ComputeRadiationForce(
    const std::array<double, 6> &_velBody,
    const std::array<double, 6> &_accelBody)
{
  // Update velocity history: newest velocity at front
  this->velHistory.push_front(_velBody);
  if (this->velHistory.size() > this->maxHistoryLen)
    this->velHistory.pop_back();

  std::array<double, 6> force{};
  force.fill(0.0);

  // Trapezoidal convolution: I_i = sum_{k=0}^{N-1} K_ij(t_k) * v_j(t-t_k) * dt
  // Trapezoidal end-point weights: w=0.5 for k=0 and k=N-1, w=1.0 otherwise
  const std::size_t N = std::min(this->velHistory.size(), this->kernel.size());
  for (std::size_t k = 0; k < N; ++k)
  {
    const double w = (k == 0 || k == N - 1) ? 0.5 : 1.0;
    const auto &Kt  = this->kernel[k];
    const auto &vel = this->velHistory[k];

    for (std::size_t row = 0; row < 6; ++row)
    {
      const std::size_t rowBase = row * 6;
      double sum = 0.0;
      for (std::size_t col = 0; col < 6; ++col)
        sum += Kt[rowBase + col] * vel[col];
      force[row] -= w * sum * this->kernelDt;
    }
  }

  // Added mass term: -A(inf) * nu_dot
  for (std::size_t row = 0; row < 6; ++row)
  {
    const std::size_t rowBase = row * 6;
    for (std::size_t col = 0; col < 6; ++col)
      force[row] -= this->addedMassInf[rowBase + col] * _accelBody[col];
  }

  return force;
}

}  // namespace usv_hydro
