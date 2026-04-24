#ifndef USV_HYDRO_HYDRO_TYPES_HH_
#define USV_HYDRO_HYDRO_TYPES_HH_

#include <array>
#include <deque>
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

  bool useLinearSeakeepingModel{false};
  std::string seakeepingCoeffsFile{};
  double seakeepingExcitationOmega{0.0};
  double seakeepingExcitationScale{1.0};

  bool useCumminsRadiation{false};
  double cumminsKernelMaxT{20.0};
  // 0 = use simulation dt automatically
  double cumminsKernelDt{0.0};

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

struct SeakeepingCoefficients
{
  std::array<double, 36> addedMass{};
  std::array<double, 36> damping{};
  std::array<double, 6> excitationRe{};
  std::array<double, 6> excitationIm{};
};

class LinearSeakeepingModel
{
  public: bool LoadFromFile(
      const std::string &_coeffFile,
      std::string *_error = nullptr);

  public: bool IsLoaded() const;

  public: bool Evaluate(
      const double _omega,
      SeakeepingCoefficients *_coeffs,
      std::string *_error = nullptr) const;

  public: std::array<double, 6> ComputeBodyWrench(
      const SeakeepingCoefficients &_coeffs,
      const std::array<double, 6> &_velocityBody,
      const std::array<double, 6> &_accelerationBody,
      const double _timeSec,
      const double _omega,
      const double _excitationScale) const;

  /// Returns only the wave excitation term (no radiation).
  /// Use this when CumminsRadiationModel handles the radiation terms.
  public: std::array<double, 6> ComputeExcitationWrench(
      const SeakeepingCoefficients &_coeffs,
      const double _timeSec,
      const double _omega,
      const double _excitationScale) const;

  public: struct FrequencySample
  {
    double omega{0.0};
    SeakeepingCoefficients coeffs;
  };

  public: const std::vector<FrequencySample> &Samples() const;

  private: std::vector<FrequencySample> samples;
};

/// Time-domain Cummins radiation model.
///
/// Implements the radiation force from the Cummins equation:
///   F_rad = -A(inf) * nu_dot - integral_0^t K(t-tau) * nu(tau) dtau
///
/// K(t) is computed from frequency-domain damping coefficients B(omega)
/// via the cosine transform:
///   K_ij(t) = (2/pi) * integral_0^inf B_ij(omega) * cos(omega*t) d_omega
///
/// A(inf) is estimated via the Ogilvie relation from K(t) and A(omega_max).
class CumminsRadiationModel
{
  /// Build the retardation kernel and A(inf) from frequency-domain samples.
  /// @param _samples  Frequency samples from LinearSeakeepingModel.
  /// @param _maxT     Maximum memory window length in seconds.
  /// @param _kernelDt Kernel time-step in seconds (must be > 0).
  /// @param _error    Optional error string on failure.
  public: bool BuildFromFrequencySamples(
      const std::vector<LinearSeakeepingModel::FrequencySample> &_samples,
      const double _maxT,
      const double _kernelDt,
      std::string *_error = nullptr);

  public: bool IsReady() const;

  /// Compute the radiation force in body frame.
  /// Internally stores _velBody in the velocity history deque.
  public: std::array<double, 6> ComputeRadiationForce(
      const std::array<double, 6> &_velBody,
      const std::array<double, 6> &_accelBody);

  /// Row-major 6x6 infinite-frequency added mass A(inf).
  private: std::array<double, 36> addedMassInf{};

  /// Retardation kernel K(t_k), kernel[k][i*6+j] = K_ij(k * kernelDt).
  private: std::vector<std::array<double, 36>> kernel;

  private: double kernelDt{0.0};

  /// Circular velocity history: velHistory[0] = current, [1] = previous, ...
  private: std::deque<std::array<double, 6>> velHistory;

  private: std::size_t maxHistoryLen{0};

  private: bool ready{false};
};

}  // namespace usv_hydro

#endif  // USV_HYDRO_HYDRO_TYPES_HH_
