#include <cmath>
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

TEST(BuoyancyModelTest, ComputesSubmergenceAndForce)
{
  usv_hydro::BuoyancyModel buoyancy(1000.0, 9.81);

  EXPECT_NEAR(buoyancy.ComputeSubmergence(-0.1, 0.5), 0.0, kEps);
  EXPECT_NEAR(buoyancy.ComputeSubmergence(0.2, 0.4), 0.5, kEps);
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
  EXPECT_NEAR(result.submergence, 0.5, kEps);
  EXPECT_NEAR(result.buoyancyWorld.Z(), 10000.0, kEps);
  EXPECT_NEAR(result.dragWorld.X(), -250.0, kEps);
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

}  // namespace
