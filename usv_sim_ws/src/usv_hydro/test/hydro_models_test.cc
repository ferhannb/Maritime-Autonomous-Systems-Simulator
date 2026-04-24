#include <array>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "usv_hydro/BuoyancyModel.hh"
#include "usv_hydro/DragModel.hh"
#include "usv_hydro/EnvironmentModel.hh"
#include "usv_hydro/HydroIntegrator.hh"
#include "usv_hydro/HydroTypes.hh"

namespace
{

constexpr double kEps = 1e-9;

struct HydroAggregate
{
  gz::math::Vector3d buoyancy{0.0, 0.0, 0.0};
  gz::math::Vector3d drag{0.0, 0.0, 0.0};
  gz::math::Vector3d totalForce{0.0, 0.0, 0.0};
  gz::math::Vector3d totalMoment{0.0, 0.0, 0.0};
};

usv_hydro::HydroCellGrid MakeBoxGrid(
    const double length,
    const double width,
    const double height,
    const int cellsX,
    const int cellsY,
    const int cellsZ)
{
  usv_hydro::HydroCellGrid grid;
  const int count = cellsX * cellsY * cellsZ;
  grid.offsets.reserve(static_cast<std::size_t>(count));

  const double minX = -0.5 * length;
  const double minY = -0.5 * width;
  const double minZ = -0.5 * height;

  for (int ix = 0; ix < cellsX; ++ix)
  {
    for (int iy = 0; iy < cellsY; ++iy)
    {
      for (int iz = 0; iz < cellsZ; ++iz)
      {
        const double x = minX + (ix + 0.5) * length / cellsX;
        const double y = minY + (iy + 0.5) * width / cellsY;
        const double z = minZ + (iz + 0.5) * height / cellsZ;
        grid.offsets.emplace_back(x, y, z);
      }
    }
  }

  grid.cellVolume = (length * width * height) / count;
  grid.cellHeightApprox = height / cellsZ;
  grid.dragAreaCell.Set(
      (width * height) / count,
      (length * height) / count,
      (length * width) / count);

  return grid;
}

HydroAggregate ComputeAggregate(
    const usv_hydro::HydroCellGrid &grid,
    const usv_hydro::HydroKinematics &kinematics,
    const usv_hydro::EnvironmentModel &environment,
    const usv_hydro::BuoyancyModel &buoyancy,
    const usv_hydro::DragModel &drag,
    const usv_hydro::HydroIntegrator &integrator)
{
  HydroAggregate aggregate;

  for (const auto &offsetLocal : grid.offsets)
  {
    const auto forces = integrator.ComputeCellForces(
        offsetLocal, grid, kinematics, environment, buoyancy, drag);
    if (!forces.submerged)
      continue;

    const auto force = forces.buoyancyWorld + forces.dragWorld;
    const auto offsetWorld = kinematics.rotationWorld.RotateVector(offsetLocal);

    aggregate.buoyancy += forces.buoyancyWorld;
    aggregate.drag += forces.dragWorld;
    aggregate.totalForce += force;
    aggregate.totalMoment += offsetWorld.Cross(force);
  }

  return aggregate;
}

TEST(HydroConfigTest, ValidatesParameters)
{
  usv_hydro::HydroConfig cfg;
  std::string error;
  EXPECT_TRUE(cfg.IsValid(&error));

  cfg.fluidDensity = 0.0;
  EXPECT_FALSE(cfg.IsValid(&error));
  EXPECT_NE(error.find("fluid_density"), std::string::npos);

  cfg = usv_hydro::HydroConfig{};
  cfg.cellsX = -1;
  EXPECT_FALSE(cfg.IsValid(&error));
  EXPECT_NE(error.find("cells_x"), std::string::npos);

  cfg = usv_hydro::HydroConfig{};
  cfg.cd.X(-0.1);
  EXPECT_FALSE(cfg.IsValid(&error));
  EXPECT_NE(error.find("cd_x"), std::string::npos);
}

TEST(HydroConfigTest, LoadsYamlProfile)
{
  const std::string seakeepingPath = "/tmp/usv_seakeeping_coeffs_test.yaml";
  std::ofstream skOut(seakeepingPath, std::ios::trunc);
  ASSERT_TRUE(skOut.good());
  skOut
      << "frequencies:\n"
      << "  - omega: 0.5\n"
      << "    added_mass: [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]\n"
      << "    damping: [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]\n";
  skOut.close();

  const std::string yamlPath = "/tmp/usv_hydro_profiles_test.yaml";
  std::ofstream out(yamlPath, std::ios::trunc);
  ASSERT_TRUE(out.good());
  out
      << "defaults:\n"
      << "  fluid_density: 1025.0\n"
      << "  cells: [12, 6, 4]\n"
      << "  drag:\n"
      << "    linear_total: [24.0, 48.0, 72.0]\n"
      << "    scale_by_cell_count: true\n"
      << "  hydrostatic_stiffness:\n"
      << "    enabled: false\n"
      << "    scale: 0.25\n"
      << "  linear_seakeeping:\n"
      << "    enabled: true\n"
      << "    coeffs_file: usv_seakeeping_coeffs_test.yaml\n"
      << "    excitation_omega: 0.8\n"
      << "    excitation_scale: 1.2\n"
      << "profiles:\n"
      << "  stiff:\n"
      << "    hydrostatic_stiffness:\n"
      << "      enabled: true\n"
      << "      scale: 0.5\n"
      << "    linear_seakeeping:\n"
      << "      excitation_scale: 0.7\n";
  out.close();

  usv_hydro::HydroConfig cfg;
  std::string error;
  ASSERT_TRUE(cfg.LoadFromFileProfile(yamlPath, "stiff", &error)) << error;

  EXPECT_NEAR(cfg.fluidDensity, 1025.0, kEps);
  EXPECT_EQ(cfg.cellsX, 12);
  EXPECT_EQ(cfg.cellsY, 6);
  EXPECT_EQ(cfg.cellsZ, 4);
  EXPECT_NEAR(cfg.linearDrag.X(), 24.0, kEps);
  EXPECT_NEAR(cfg.linearDrag.Y(), 48.0, kEps);
  EXPECT_NEAR(cfg.linearDrag.Z(), 72.0, kEps);
  EXPECT_TRUE(cfg.scaleLinearDragByCellCount);
  EXPECT_TRUE(cfg.useHydrostaticStiffnessMatrix);
  EXPECT_NEAR(cfg.hydrostaticStiffnessScale, 0.5, kEps);
  EXPECT_TRUE(cfg.useLinearSeakeepingModel);
  EXPECT_EQ(cfg.seakeepingCoeffsFile, seakeepingPath);
  EXPECT_NEAR(cfg.seakeepingExcitationOmega, 0.8, kEps);
  EXPECT_NEAR(cfg.seakeepingExcitationScale, 0.7, kEps);
}

TEST(LinearSeakeepingModelTest, InterpolatesAndComputesBodyWrench)
{
  const std::string yamlPath = "/tmp/usv_linear_seakeeping_test.yaml";
  std::ofstream out(yamlPath, std::ios::trunc);
  ASSERT_TRUE(out.good());
  out
      << "frequencies:\n"
      << "  - omega: 0.5\n"
      << "    added_mass: [10, 0, 0, 0, 0, 0, 0, 11, 0, 0, 0, 0, 0, 0, 12, 0, 0, 0, 0, 0, 0, 13, 0, 0, 0, 0, 0, 0, 14, 0, 0, 0, 0, 0, 0, 15]\n"
      << "    damping: [1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 6]\n"
      << "    excitation_re: [100, 0, 0, 0, 0, 0]\n"
      << "    excitation_im: [0, 0, 0, 0, 0, 0]\n"
      << "  - omega: 1.5\n"
      << "    added_mass: [20, 0, 0, 0, 0, 0, 0, 21, 0, 0, 0, 0, 0, 0, 22, 0, 0, 0, 0, 0, 0, 23, 0, 0, 0, 0, 0, 0, 24, 0, 0, 0, 0, 0, 0, 25]\n"
      << "    damping: [2, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 7]\n"
      << "    excitation_re: [200, 0, 0, 0, 0, 0]\n"
      << "    excitation_im: [0, 0, 0, 0, 0, 0]\n";
  out.close();

  usv_hydro::LinearSeakeepingModel model;
  std::string error;
  ASSERT_TRUE(model.LoadFromFile(yamlPath, &error)) << error;
  ASSERT_TRUE(model.IsLoaded());

  usv_hydro::SeakeepingCoefficients coeffs;
  ASSERT_TRUE(model.Evaluate(1.0, &coeffs, &error)) << error;

  EXPECT_NEAR(coeffs.addedMass[0], 15.0, kEps);
  EXPECT_NEAR(coeffs.addedMass[7], 16.0, kEps);
  EXPECT_NEAR(coeffs.damping[0], 1.5, kEps);
  EXPECT_NEAR(coeffs.excitationRe[0], 150.0, kEps);

  const std::array<double, 6> velocity{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const std::array<double, 6> acceleration{0.2, 0.0, 0.0, 0.0, 0.0, 0.0};
  const auto wrench = model.ComputeBodyWrench(
      coeffs, velocity, acceleration, 0.0, 1.0, 1.0);

  // -A*u_dot - B*u + Re{F_exc*exp(iwt)} at t=0.
  EXPECT_NEAR(wrench[0], -(15.0 * 0.2 + 1.5 * 1.0) + 150.0, kEps);
  EXPECT_NEAR(wrench[1], 0.0, kEps);
}

TEST(BuoyancyModelTest, ComputesSubmergenceAndForce)
{
  usv_hydro::BuoyancyModel buoyancy(1000.0, 9.81);

  EXPECT_NEAR(buoyancy.ComputeSubmergence(-0.3, 0.5), 0.0, kEps);
  EXPECT_NEAR(buoyancy.ComputeSubmergence(-0.2, 0.4), 0.0, kEps);
  EXPECT_NEAR(buoyancy.ComputeSubmergence(0.0, 0.4), 0.5, kEps);
  EXPECT_NEAR(buoyancy.ComputeSubmergence(0.2, 0.4), 1.0, kEps);
  EXPECT_NEAR(buoyancy.ComputeSubmergence(10.0, 0.4), 1.0, kEps);

  const auto force = buoyancy.ComputeForceWorld(2.0, 0.25);
  EXPECT_NEAR(force.X(), 0.0, kEps);
  EXPECT_NEAR(force.Y(), 0.0, kEps);
  EXPECT_NEAR(force.Z(), 1000.0 * 9.81 * 0.5, kEps);
}

TEST(DragModelTest, ComputesLinearAndQuadraticDrag)
{
  const usv_hydro::DragModel drag(
      1000.0,
      gz::math::Vector3d(1.0, 2.0, 0.5),
      gz::math::Vector3d(10.0, 20.0, 30.0));

  const gz::math::Vector3d velBody(2.0, -3.0, 0.5);
  const gz::math::Vector3d area(0.2, 0.3, 0.4);

  const auto force = drag.ComputeForceBody(velBody, area);

  const double fx = -10.0 * 2.0 - 0.5 * 1000.0 * 1.0 * 0.2 * std::abs(2.0) * 2.0;
  const double fy =
      -20.0 * -3.0 - 0.5 * 1000.0 * 2.0 * 0.3 * std::abs(-3.0) * -3.0;
  const double fz =
      -30.0 * 0.5 - 0.5 * 1000.0 * 0.5 * 0.4 * std::abs(0.5) * 0.5;

  EXPECT_NEAR(force.X(), fx, kEps);
  EXPECT_NEAR(force.Y(), fy, kEps);
  EXPECT_NEAR(force.Z(), fz, kEps);
}

TEST(EnvironmentModelTest, ComputesDepthAndRelativeVelocity)
{
  usv_hydro::HydroConfig cfg;
  cfg.waterLevel = 1.5;
  cfg.fluidDensity = 997.0;
  cfg.gravity = 9.7;
  cfg.currentVelocity = gz::math::Vector3d(0.5, -0.2, 0.0);

  const usv_hydro::EnvironmentModel env(cfg);

  EXPECT_NEAR(env.DepthAt(0.5), 1.0, kEps);
  EXPECT_NEAR(env.FluidDensity(), 997.0, kEps);
  EXPECT_NEAR(env.Gravity(), 9.7, kEps);

  const auto rel = env.RelativeVelocityWorld(gz::math::Vector3d(1.0, 0.0, 0.0));
  EXPECT_NEAR(rel.X(), 0.5, kEps);
  EXPECT_NEAR(rel.Y(), 0.2, kEps);
  EXPECT_NEAR(rel.Z(), 0.0, kEps);
}

TEST(HydroIntegratorTest, ComputesCellForcesWhenSubmerged)
{
  usv_hydro::HydroCellGrid grid;
  grid.cellVolume = 2.0;
  grid.cellHeightApprox = 1.0;
  grid.dragAreaCell = gz::math::Vector3d(1.0, 1.0, 1.0);

  usv_hydro::HydroKinematics kinematics;
  kinematics.positionWorld = gz::math::Vector3d(0.0, 0.0, 0.0);
  kinematics.rotationWorld = gz::math::Quaterniond(0.0, 0.0, 0.0);
  kinematics.linearVelocityWorld = gz::math::Vector3d(1.0, 0.0, 0.0);
  kinematics.angularVelocityWorld = gz::math::Vector3d(0.0, 0.0, 0.0);

  usv_hydro::HydroConfig cfg;
  cfg.waterLevel = 0.5;
  cfg.fluidDensity = 1000.0;
  cfg.gravity = 10.0;
  cfg.currentVelocity = gz::math::Vector3d(0.0, 0.0, 0.0);
  cfg.cd = gz::math::Vector3d(1.0, 1.0, 1.0);
  cfg.linearDrag = gz::math::Vector3d(0.0, 0.0, 0.0);

  const usv_hydro::EnvironmentModel env(cfg);
  const usv_hydro::BuoyancyModel buoyancy(cfg.fluidDensity, cfg.gravity);
  const usv_hydro::DragModel drag(cfg.fluidDensity, cfg.cd, cfg.linearDrag);
  const usv_hydro::HydroIntegrator integrator;

  const auto result = integrator.ComputeCellForces(
      gz::math::Vector3d(0.0, 0.0, 0.0), grid, kinematics, env, buoyancy, drag);

  EXPECT_TRUE(result.submerged);
  EXPECT_NEAR(result.submergence, 1.0, kEps);
  EXPECT_NEAR(result.buoyancyWorld.Z(), 20000.0, kEps);
  EXPECT_NEAR(result.dragWorld.X(), -500.0, kEps);
  EXPECT_NEAR(result.dragWorld.Y(), 0.0, kEps);
  EXPECT_NEAR(result.dragWorld.Z(), 0.0, kEps);
}

TEST(HydroIntegratorTest, ReturnsZeroForcesOutOfWater)
{
  usv_hydro::HydroCellGrid grid;
  grid.cellVolume = 2.0;
  grid.cellHeightApprox = 1.0;
  grid.dragAreaCell = gz::math::Vector3d(1.0, 1.0, 1.0);

  usv_hydro::HydroKinematics kinematics;
  kinematics.positionWorld = gz::math::Vector3d(0.0, 0.0, 2.0);
  kinematics.rotationWorld = gz::math::Quaterniond(0.0, 0.0, 0.0);

  usv_hydro::HydroConfig cfg;
  cfg.waterLevel = 0.5;
  const usv_hydro::EnvironmentModel env(cfg);
  const usv_hydro::BuoyancyModel buoyancy(cfg.fluidDensity, cfg.gravity);
  const usv_hydro::DragModel drag(cfg.fluidDensity, cfg.cd, cfg.linearDrag);
  const usv_hydro::HydroIntegrator integrator;

  const auto result = integrator.ComputeCellForces(
      gz::math::Vector3d(0.0, 0.0, 0.0), grid, kinematics, env, buoyancy, drag);

  EXPECT_FALSE(result.submerged);
  EXPECT_NEAR(result.submergence, 0.0, kEps);
  EXPECT_NEAR(result.buoyancyWorld.Length(), 0.0, kEps);
  EXPECT_NEAR(result.dragWorld.Length(), 0.0, kEps);
}

TEST(HydroScenarioTest, HydrostaticEquilibriumSettlesNearAnalyticalDraft)
{
  constexpr double hullLength = 2.2;
  constexpr double hullWidth = 1.1;
  constexpr double hullHeight = 0.5;
  constexpr double mass = 350.0;
  constexpr double rho = 1000.0;
  constexpr double gravity = 9.81;
  constexpr double dt = 0.01;
  constexpr int steps = 4000;

  const auto grid = MakeBoxGrid(hullLength, hullWidth, hullHeight, 14, 8, 4);

  usv_hydro::HydroConfig cfg;
  cfg.waterLevel = 0.0;
  cfg.fluidDensity = rho;
  cfg.gravity = gravity;
  cfg.cd = gz::math::Vector3d(0.9, 1.3, 2.0);
  cfg.linearDrag = gz::math::Vector3d(18.0, 55.0, 80.0);

  const usv_hydro::EnvironmentModel env(cfg);
  const usv_hydro::BuoyancyModel buoyancy(rho, gravity);
  const usv_hydro::DragModel drag(rho, cfg.cd, cfg.linearDrag);
  const usv_hydro::HydroIntegrator integrator;

  usv_hydro::HydroKinematics kinematics;
  kinematics.rotationWorld = gz::math::Quaterniond(0.0, 0.0, 0.0);

  double z = 0.40;
  double vz = 0.0;

  const auto netVerticalForce = [&](const double zEval) {
    usv_hydro::HydroKinematics kinematicsAtZ;
    kinematicsAtZ.positionWorld = gz::math::Vector3d(0.0, 0.0, zEval);
    kinematicsAtZ.rotationWorld = gz::math::Quaterniond(0.0, 0.0, 0.0);
    kinematicsAtZ.linearVelocityWorld = gz::math::Vector3d(0.0, 0.0, 0.0);
    kinematicsAtZ.angularVelocityWorld = gz::math::Vector3d(0.0, 0.0, 0.0);
    const auto aggregate =
        ComputeAggregate(grid, kinematicsAtZ, env, buoyancy, drag, integrator);
    return aggregate.totalForce.Z() - mass * gravity;
  };

  double zLow = -0.5;
  double zHigh = 0.5;
  ASSERT_GT(netVerticalForce(zLow), 0.0);
  ASSERT_LT(netVerticalForce(zHigh), 0.0);
  for (int i = 0; i < 80; ++i)
  {
    const double zMid = 0.5 * (zLow + zHigh);
    if (netVerticalForce(zMid) > 0.0)
      zLow = zMid;
    else
      zHigh = zMid;
  }
  const double discretizedEquilibriumZ = 0.5 * (zLow + zHigh);

  for (int i = 0; i < steps; ++i)
  {
    kinematics.positionWorld = gz::math::Vector3d(0.0, 0.0, z);
    kinematics.linearVelocityWorld = gz::math::Vector3d(0.0, 0.0, vz);
    kinematics.angularVelocityWorld = gz::math::Vector3d(0.0, 0.0, 0.0);

    const auto aggregate =
        ComputeAggregate(grid, kinematics, env, buoyancy, drag, integrator);
    const double netFz = aggregate.totalForce.Z() - mass * gravity;
    const double az = netFz / mass;
    vz += az * dt;
    z += vz * dt;
  }

  const double analyticalDraft = mass / (rho * hullLength * hullWidth);
  const double expectedCenterZ = 0.5 * hullHeight - analyticalDraft;
  EXPECT_NEAR(
      discretizedEquilibriumZ,
      expectedCenterZ,
      0.5 * grid.cellHeightApprox + 0.02);
  EXPECT_NEAR(z, discretizedEquilibriumZ, 0.02);
}

TEST(HydroScenarioTest, SurgeVelocityDecaysWithoutThrust)
{
  constexpr double mass = 350.0;
  constexpr double dt = 0.01;
  constexpr int steps = 2000;

  const auto grid = MakeBoxGrid(2.2, 1.1, 0.5, 14, 8, 4);

  usv_hydro::HydroConfig cfg;
  cfg.waterLevel = 0.0;
  cfg.fluidDensity = 1000.0;
  cfg.gravity = 9.81;
  cfg.cd = gz::math::Vector3d(0.9, 1.3, 2.0);
  cfg.linearDrag = gz::math::Vector3d(18.0, 55.0, 80.0);

  const usv_hydro::EnvironmentModel env(cfg);
  const usv_hydro::BuoyancyModel buoyancy(cfg.fluidDensity, cfg.gravity);
  const usv_hydro::DragModel drag(cfg.fluidDensity, cfg.cd, cfg.linearDrag);
  const usv_hydro::HydroIntegrator integrator;

  usv_hydro::HydroKinematics kinematics;
  kinematics.positionWorld = gz::math::Vector3d(0.0, 0.0, 0.1053719);
  kinematics.rotationWorld = gz::math::Quaterniond(0.0, 0.0, 0.0);
  kinematics.angularVelocityWorld = gz::math::Vector3d(0.0, 0.0, 0.0);

  double vx = 2.0;
  for (int i = 0; i < steps; ++i)
  {
    kinematics.linearVelocityWorld = gz::math::Vector3d(vx, 0.0, 0.0);
    const auto aggregate =
        ComputeAggregate(grid, kinematics, env, buoyancy, drag, integrator);

    const double ax = aggregate.totalForce.X() / mass;
    const double vxNext = vx + ax * dt;
    EXPECT_LE(vxNext, vx + 1e-6);
    vx = vxNext;
  }

  EXPECT_LT(vx, 0.2);
  EXPECT_GT(vx, 0.0);
}

TEST(HydroScenarioTest, DragCancelsAtCurrentMatchedVelocity)
{
  const auto grid = MakeBoxGrid(2.2, 1.1, 0.5, 14, 8, 4);

  usv_hydro::HydroConfig cfg;
  cfg.waterLevel = 0.0;
  cfg.fluidDensity = 1000.0;
  cfg.gravity = 9.81;
  cfg.currentVelocity = gz::math::Vector3d(1.2, -0.4, 0.0);

  const usv_hydro::EnvironmentModel env(cfg);
  const usv_hydro::BuoyancyModel buoyancy(cfg.fluidDensity, cfg.gravity);
  const usv_hydro::DragModel drag(cfg.fluidDensity, cfg.cd, cfg.linearDrag);
  const usv_hydro::HydroIntegrator integrator;

  usv_hydro::HydroKinematics kinematics;
  kinematics.positionWorld = gz::math::Vector3d(0.0, 0.0, 0.1053719);
  kinematics.rotationWorld = gz::math::Quaterniond(0.0, 0.0, 0.0);
  kinematics.linearVelocityWorld = cfg.currentVelocity;
  kinematics.angularVelocityWorld = gz::math::Vector3d(0.0, 0.0, 0.0);

  const auto aggregate =
      ComputeAggregate(grid, kinematics, env, buoyancy, drag, integrator);

  EXPECT_NEAR(aggregate.drag.X(), 0.0, 1e-9);
  EXPECT_NEAR(aggregate.drag.Y(), 0.0, 1e-9);
  EXPECT_NEAR(aggregate.drag.Z(), 0.0, 1e-9);
}

TEST(HydroScenarioTest, SymmetryKeepsRollPitchMomentsNearZero)
{
  const auto grid = MakeBoxGrid(2.2, 1.1, 0.5, 14, 8, 4);

  usv_hydro::HydroConfig cfg;
  cfg.waterLevel = 0.0;
  cfg.fluidDensity = 1000.0;
  cfg.gravity = 9.81;

  const usv_hydro::EnvironmentModel env(cfg);
  const usv_hydro::BuoyancyModel buoyancy(cfg.fluidDensity, cfg.gravity);
  const usv_hydro::DragModel drag(cfg.fluidDensity, cfg.cd, cfg.linearDrag);
  const usv_hydro::HydroIntegrator integrator;

  usv_hydro::HydroKinematics kinematics;
  kinematics.positionWorld = gz::math::Vector3d(0.0, 0.0, 0.1053719);
  kinematics.rotationWorld = gz::math::Quaterniond(0.0, 0.0, 0.0);
  kinematics.linearVelocityWorld = gz::math::Vector3d(0.0, 0.0, 0.0);
  kinematics.angularVelocityWorld = gz::math::Vector3d(0.0, 0.0, 0.0);

  const auto aggregate =
      ComputeAggregate(grid, kinematics, env, buoyancy, drag, integrator);

  EXPECT_NEAR(aggregate.totalMoment.X(), 0.0, 1e-6);
  EXPECT_NEAR(aggregate.totalMoment.Y(), 0.0, 1e-6);
}

TEST(HydroScenarioTest, HydrostaticStiffnessMatchesBoxTheory)
{
  constexpr double length = 2.2;
  constexpr double width = 1.1;
  constexpr double height = 0.5;
  constexpr double mass = 350.0;
  constexpr double rho = 1000.0;
  constexpr double gravity = 9.81;

  const auto grid = MakeBoxGrid(length, width, height, 28, 16, 10);

  usv_hydro::HydroConfig cfg;
  cfg.waterLevel = 0.0;
  cfg.fluidDensity = rho;
  cfg.gravity = gravity;
  cfg.currentVelocity = gz::math::Vector3d(0.0, 0.0, 0.0);

  const usv_hydro::EnvironmentModel env(cfg);
  const usv_hydro::BuoyancyModel buoyancy(cfg.fluidDensity, cfg.gravity);
  const usv_hydro::DragModel drag(cfg.fluidDensity, cfg.cd, cfg.linearDrag);
  const usv_hydro::HydroIntegrator integrator;

  const double draft = mass / (rho * length * width);
  const double zEq = 0.5 * height - draft;

  const auto aggregateAt = [&](const double z, const double roll, const double pitch)
  {
    usv_hydro::HydroKinematics kinematics;
    kinematics.positionWorld = gz::math::Vector3d(0.0, 0.0, z);
    kinematics.rotationWorld = gz::math::Quaterniond(roll, pitch, 0.0);
    kinematics.linearVelocityWorld = gz::math::Vector3d(0.0, 0.0, 0.0);
    kinematics.angularVelocityWorld = gz::math::Vector3d(0.0, 0.0, 0.0);
    return ComputeAggregate(grid, kinematics, env, buoyancy, drag, integrator);
  };

  const double dz = 0.01;
  const double dTheta = 1e-3;

  const double fzPlus = aggregateAt(zEq + dz, 0.0, 0.0).totalForce.Z();
  const double fzMinus = aggregateAt(zEq - dz, 0.0, 0.0).totalForce.Z();
  const double k33Num = -(fzPlus - fzMinus) / (2.0 * dz);

  const double fzRollPlus = aggregateAt(zEq, dTheta, 0.0).totalForce.Z();
  const double fzRollMinus = aggregateAt(zEq, -dTheta, 0.0).totalForce.Z();
  const double k34Num = -(fzRollPlus - fzRollMinus) / (2.0 * dTheta);

  const double fzPitchPlus = aggregateAt(zEq, 0.0, dTheta).totalForce.Z();
  const double fzPitchMinus = aggregateAt(zEq, 0.0, -dTheta).totalForce.Z();
  const double k35Num = -(fzPitchPlus - fzPitchMinus) / (2.0 * dTheta);

  const double mxPlus = aggregateAt(zEq, dTheta, 0.0).totalMoment.X();
  const double mxMinus = aggregateAt(zEq, -dTheta, 0.0).totalMoment.X();
  const double k44Num = -(mxPlus - mxMinus) / (2.0 * dTheta);

  const double mxHeavePlus = aggregateAt(zEq + dz, 0.0, 0.0).totalMoment.X();
  const double mxHeaveMinus = aggregateAt(zEq - dz, 0.0, 0.0).totalMoment.X();
  const double k43Num = -(mxHeavePlus - mxHeaveMinus) / (2.0 * dz);

  const double mxPitchPlus = aggregateAt(zEq, 0.0, dTheta).totalMoment.X();
  const double mxPitchMinus = aggregateAt(zEq, 0.0, -dTheta).totalMoment.X();
  const double k45Num = -(mxPitchPlus - mxPitchMinus) / (2.0 * dTheta);

  const double myPlus = aggregateAt(zEq, 0.0, dTheta).totalMoment.Y();
  const double myMinus = aggregateAt(zEq, 0.0, -dTheta).totalMoment.Y();
  const double k55Num = -(myPlus - myMinus) / (2.0 * dTheta);

  const double myHeavePlus = aggregateAt(zEq + dz, 0.0, 0.0).totalMoment.Y();
  const double myHeaveMinus = aggregateAt(zEq - dz, 0.0, 0.0).totalMoment.Y();
  const double k53Num = -(myHeavePlus - myHeaveMinus) / (2.0 * dz);

  const double myRollPlus = aggregateAt(zEq, dTheta, 0.0).totalMoment.Y();
  const double myRollMinus = aggregateAt(zEq, -dTheta, 0.0).totalMoment.Y();
  const double k54Num = -(myRollPlus - myRollMinus) / (2.0 * dTheta);

  const double displacedVolume = mass / rho;
  const double waterplaneArea = length * width;
  const double iXX = length * std::pow(width, 3) / 12.0;
  const double iYY = width * std::pow(length, 3) / 12.0;
  const double zBRelCG = 0.5 * (draft - height);
  const double gmT = iXX / displacedVolume + zBRelCG;
  const double gmL = iYY / displacedVolume + zBRelCG;

  const double k33Theory = rho * gravity * waterplaneArea;
  const double k44Theory = rho * gravity * displacedVolume * gmT;
  const double k55Theory = rho * gravity * displacedVolume * gmL;

  EXPECT_GT(k33Num, 0.0);
  EXPECT_GT(k44Num, 0.0);
  EXPECT_GT(k55Num, 0.0);

  EXPECT_NEAR(k33Num, k33Theory, 0.10 * k33Theory);
  EXPECT_NEAR(k44Num, k44Theory, 0.15 * k44Theory);
  EXPECT_NEAR(k55Num, k55Theory, 0.15 * k55Theory);

  // For the symmetric box case at centered COG, cross-couplings should vanish.
  EXPECT_NEAR(k34Num, 0.0, 200.0);
  EXPECT_NEAR(k35Num, 0.0, 200.0);
  EXPECT_NEAR(k43Num, 0.0, 200.0);
  EXPECT_NEAR(k53Num, 0.0, 200.0);
  EXPECT_NEAR(k45Num, 0.0, 200.0);
  EXPECT_NEAR(k54Num, 0.0, 200.0);

  // Hydrostatic matrix symmetry check.
  EXPECT_NEAR(k34Num, k43Num, 100.0);
  EXPECT_NEAR(k35Num, k53Num, 100.0);
  EXPECT_NEAR(k45Num, k54Num, 100.0);
}

// ---------------------------------------------------------------------------
// Helpers shared by new tests
// ---------------------------------------------------------------------------

static std::vector<usv_hydro::LinearSeakeepingModel::FrequencySample>
MakeDiagonalSamples(const std::vector<double> &omegas,
                    const double bDiag,
                    const double aDiag)
{
  std::vector<usv_hydro::LinearSeakeepingModel::FrequencySample> samples;
  samples.reserve(omegas.size());
  for (const double om : omegas)
  {
    usv_hydro::LinearSeakeepingModel::FrequencySample s;
    s.omega = om;
    for (int i = 0; i < 6; ++i)
    {
      s.coeffs.damping[i * 6 + i] = bDiag;
      s.coeffs.addedMass[i * 6 + i] = aDiag;
    }
    samples.push_back(s);
  }
  return samples;
}

// ---------------------------------------------------------------------------
// CumminsRadiationModel tests
// ---------------------------------------------------------------------------

TEST(CumminsRadiationModelTest, BuildFromTwoFrequencySamples)
{
  const auto samples = MakeDiagonalSamples({0.5, 1.5}, 50.0, 100.0);

  usv_hydro::CumminsRadiationModel model;
  std::string error;
  ASSERT_TRUE(model.BuildFromFrequencySamples(samples, 10.0, 0.1, &error))
      << error;
  EXPECT_TRUE(model.IsReady());
}

TEST(CumminsRadiationModelTest, RequiresAtLeastTwoSamples)
{
  const auto samples = MakeDiagonalSamples({1.0}, 50.0, 100.0);

  usv_hydro::CumminsRadiationModel model;
  std::string error;
  EXPECT_FALSE(model.BuildFromFrequencySamples(samples, 10.0, 0.1, &error));
  EXPECT_FALSE(model.IsReady());
}

TEST(CumminsRadiationModelTest, ZeroVelocityYieldsZeroMemoryForce)
{
  const auto samples = MakeDiagonalSamples({0.5, 1.5}, 50.0, 100.0);

  usv_hydro::CumminsRadiationModel model;
  std::string error;
  ASSERT_TRUE(model.BuildFromFrequencySamples(samples, 10.0, 0.1, &error))
      << error;

  const std::array<double, 6> zeroVel{};
  const std::array<double, 6> zeroAccel{};

  // Fill history with zeros, then verify force remains zero
  for (int i = 0; i < 120; ++i)
  {
    const auto force = model.ComputeRadiationForce(zeroVel, zeroAccel);
    for (int j = 0; j < 6; ++j)
      EXPECT_NEAR(force[j], 0.0, kEps) << "step " << i << " dof " << j;
  }
}

TEST(CumminsRadiationModelTest, AddedMassTermNonZeroForNonZeroAccel)
{
  const auto samples = MakeDiagonalSamples({0.5, 1.5}, 50.0, 150.0);

  usv_hydro::CumminsRadiationModel model;
  std::string error;
  ASSERT_TRUE(model.BuildFromFrequencySamples(samples, 10.0, 0.1, &error))
      << error;

  const std::array<double, 6> zeroVel{};
  const std::array<double, 6> accel{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  const auto force = model.ComputeRadiationForce(zeroVel, accel);
  // F[0] ≈ -A_inf_00 * accel[0]; A_inf_00 > 0 → force[0] < 0
  EXPECT_LT(force[0], 0.0);
}

TEST(CumminsRadiationModelTest, MemoryForceBuildsWithNonZeroVelocityHistory)
{
  // Use 10 omega samples for a broader-band kernel
  const std::vector<double> omegas = {0.2, 0.4, 0.6, 0.8, 1.0,
                                       1.2, 1.4, 1.6, 1.8, 2.0};
  const auto samples = MakeDiagonalSamples(omegas, 100.0, 150.0);

  usv_hydro::CumminsRadiationModel model;
  std::string error;
  ASSERT_TRUE(model.BuildFromFrequencySamples(samples, 20.0, 0.05, &error))
      << error;

  const std::array<double, 6> vel{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const std::array<double, 6> zeroAccel{};

  // After filling history with constant unit surge velocity, radiation force
  // in surge (index 0) should be non-trivially negative (radiation damping)
  std::array<double, 6> lastForce{};
  for (int i = 0; i < 50; ++i)
    lastForce = model.ComputeRadiationForce(vel, zeroAccel);

  EXPECT_LT(lastForce[0], 0.0);
}

TEST(CumminsRadiationModelTest, IsValidatedByConfig)
{
  const std::string skPath = "/tmp/usv_cummins_sk_test.yaml";
  std::ofstream out(skPath, std::ios::trunc);
  ASSERT_TRUE(out.good());
  out << "frequencies:\n"
      << "  - omega: 0.5\n"
      << "    added_mass: [1,0,0,0,0,0, 0,1,0,0,0,0, 0,0,1,0,0,0,"
         " 0,0,0,1,0,0, 0,0,0,0,1,0, 0,0,0,0,0,1]\n"
      << "    damping: [1,0,0,0,0,0, 0,1,0,0,0,0, 0,0,1,0,0,0,"
         " 0,0,0,1,0,0, 0,0,0,0,1,0, 0,0,0,0,0,1]\n"
      << "  - omega: 1.5\n"
      << "    added_mass: [2,0,0,0,0,0, 0,2,0,0,0,0, 0,0,2,0,0,0,"
         " 0,0,0,2,0,0, 0,0,0,0,2,0, 0,0,0,0,0,2]\n"
      << "    damping: [2,0,0,0,0,0, 0,2,0,0,0,0, 0,0,2,0,0,0,"
         " 0,0,0,2,0,0, 0,0,0,0,2,0, 0,0,0,0,0,2]\n";
  out.close();

  usv_hydro::HydroConfig cfg;
  cfg.useLinearSeakeepingModel = true;
  cfg.seakeepingCoeffsFile = skPath;
  cfg.useCumminsRadiation = true;
  cfg.cumminsKernelMaxT = 15.0;
  cfg.cumminsKernelDt = 0.0;

  std::string error;
  EXPECT_TRUE(cfg.IsValid(&error)) << error;

  // Cummins without seakeeping should fail validation
  usv_hydro::HydroConfig bad;
  bad.useCumminsRadiation = true;
  bad.cumminsKernelMaxT = 15.0;
  EXPECT_FALSE(bad.IsValid(&error));
  EXPECT_NE(error.find("cummins"), std::string::npos);

  // Cummins with negative max_t should fail
  usv_hydro::HydroConfig bad2;
  bad2.useLinearSeakeepingModel = true;
  bad2.seakeepingCoeffsFile = skPath;
  bad2.useCumminsRadiation = true;
  bad2.cumminsKernelMaxT = -1.0;
  EXPECT_FALSE(bad2.IsValid(&error));
}

// ---------------------------------------------------------------------------
// 6x6 hydrostatic stiffness tests
// ---------------------------------------------------------------------------

TEST(HydroScenarioTest, StiffnessSurgeSwayNearZeroForSymmetricHull)
{
  // For a box hull with flat bottom, a horizontal shift does not change the
  // displaced volume → surge/sway restoring stiffness should be ≈ 0.
  constexpr double length = 2.2;
  constexpr double width  = 1.1;
  constexpr double height = 0.5;
  constexpr double mass   = 350.0;
  constexpr double rho    = 1000.0;
  constexpr double gravity = 9.81;

  const auto grid = MakeBoxGrid(length, width, height, 28, 16, 10);

  usv_hydro::HydroConfig cfg;
  cfg.waterLevel   = 0.0;
  cfg.fluidDensity = rho;
  cfg.gravity      = gravity;

  const usv_hydro::EnvironmentModel env(cfg);
  const usv_hydro::BuoyancyModel buoyancy(rho, gravity);
  const usv_hydro::DragModel drag(rho, cfg.cd, cfg.linearDrag);
  const usv_hydro::HydroIntegrator integrator;

  const double draft = mass / (rho * length * width);
  const double zEq   = 0.5 * height - draft;
  const double ds    = 0.01;

  const auto aggAt = [&](double x, double y, double z)
  {
    usv_hydro::HydroKinematics kin;
    kin.positionWorld = gz::math::Vector3d(x, y, z);
    kin.rotationWorld = gz::math::Quaterniond(0.0, 0.0, 0.0);
    kin.linearVelocityWorld  = gz::math::Vector3d(0.0, 0.0, 0.0);
    kin.angularVelocityWorld = gz::math::Vector3d(0.0, 0.0, 0.0);
    return ComputeAggregate(grid, kin, env, buoyancy, drag, integrator);
  };

  // C11 (surge column, surge row): dFx/dx_surge
  const double fxP = aggAt( ds, 0, zEq).totalForce.X();
  const double fxM = aggAt(-ds, 0, zEq).totalForce.X();
  const double c11 = -(fxP - fxM) / (2.0 * ds);
  EXPECT_NEAR(c11, 0.0, 5.0);

  // C22 (sway column, sway row): dFy/dy_sway
  const double fyP = aggAt(0,  ds, zEq).totalForce.Y();
  const double fyM = aggAt(0, -ds, zEq).totalForce.Y();
  const double c22 = -(fyP - fyM) / (2.0 * ds);
  EXPECT_NEAR(c22, 0.0, 5.0);
}

TEST(HydroScenarioTest, Stiffness6x6Reciprocity)
{
  // Verify C35 ≈ C53 (heave-pitch / pitch-heave coupling reciprocity)
  // and C34 ≈ C43 for the symmetric box hull.
  constexpr double length = 2.2;
  constexpr double width  = 1.1;
  constexpr double height = 0.5;
  constexpr double mass   = 350.0;
  constexpr double rho    = 1000.0;
  constexpr double gravity = 9.81;

  const auto grid = MakeBoxGrid(length, width, height, 28, 16, 10);

  usv_hydro::HydroConfig cfg;
  cfg.waterLevel   = 0.0;
  cfg.fluidDensity = rho;
  cfg.gravity      = gravity;

  const usv_hydro::EnvironmentModel env(cfg);
  const usv_hydro::BuoyancyModel buoyancy(rho, gravity);
  const usv_hydro::DragModel drag(rho, cfg.cd, cfg.linearDrag);
  const usv_hydro::HydroIntegrator integrator;

  const double draft  = mass / (rho * length * width);
  const double zEq    = 0.5 * height - draft;
  const double dz     = 0.01;
  const double dTheta = 1e-3;

  const auto aggAt = [&](double z, double roll, double pitch)
  {
    usv_hydro::HydroKinematics kin;
    kin.positionWorld = gz::math::Vector3d(0.0, 0.0, z);
    kin.rotationWorld = gz::math::Quaterniond(roll, pitch, 0.0);
    kin.linearVelocityWorld  = gz::math::Vector3d(0.0, 0.0, 0.0);
    kin.angularVelocityWorld = gz::math::Vector3d(0.0, 0.0, 0.0);
    return ComputeAggregate(grid, kin, env, buoyancy, drag, integrator);
  };

  // C35 = -dFz/d(pitch)  (heave force per unit pitch)
  const double fzPP = aggAt(zEq, 0.0,  dTheta).totalForce.Z();
  const double fzPM = aggAt(zEq, 0.0, -dTheta).totalForce.Z();
  const double c35  = -(fzPP - fzPM) / (2.0 * dTheta);

  // C53 = -dMy/d(heave)  (pitch moment per unit heave)
  const double myHP = aggAt(zEq + dz, 0.0, 0.0).totalMoment.Y();
  const double myHM = aggAt(zEq - dz, 0.0, 0.0).totalMoment.Y();
  const double c53  = -(myHP - myHM) / (2.0 * dz);

  EXPECT_NEAR(c35, c53, 150.0);

  // C34 = -dFz/d(roll)  (heave force per unit roll) — near zero for symmetric hull
  const double fzRP = aggAt(zEq,  dTheta, 0.0).totalForce.Z();
  const double fzRM = aggAt(zEq, -dTheta, 0.0).totalForce.Z();
  const double c34  = -(fzRP - fzRM) / (2.0 * dTheta);

  // C43 = -dMx/d(heave)  (roll moment per unit heave) — near zero for symmetric hull
  const double mxHP = aggAt(zEq + dz, 0.0, 0.0).totalMoment.X();
  const double mxHM = aggAt(zEq - dz, 0.0, 0.0).totalMoment.X();
  const double c43  = -(mxHP - mxHM) / (2.0 * dz);

  EXPECT_NEAR(c34, c43, 100.0);
}

}  // namespace
