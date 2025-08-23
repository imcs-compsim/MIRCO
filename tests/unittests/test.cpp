#include <gtest/gtest.h>
#include <stdlib.h>

#include <vector>

#include "../../src/mirco_inputparameters_kokkos.h"
#include "../../src/mirco_kokkostypes_kokkos.h"
#include "../../src/mirco_nonlinearsolver_kokkos.h"
#include "../../src/mirco_topology_kokkos.h"
#include "../../src/mirco_utils_kokkos.h"
#include "../../src/mirco_warmstart_kokkos.h"

TEST(NonlinearSolverTest, primalvariable)
{
  MIRCO::ViewMatrix_h matrix_h("matrix_h", 9, 9);
  matrix_h(0, 0) = 0.00381971863420549;
  matrix_h(1, 0) = 0.0020;
  matrix_h(2, 0) = 0.000965167479061995;
  matrix_h(3, 0) = 0.0020;
  matrix_h(4, 0) = 0.00138032073697570;
  matrix_h(5, 0) = 0.000861397758772237;
  matrix_h(6, 0) = 0.000965167479061995;
  matrix_h(7, 0) = 0.000861397758772237;
  matrix_h(8, 0) = 0.000678804493543927;

  matrix_h(0, 1) = 0.002;
  matrix_h(1, 1) = 0.00381971863420549;
  matrix_h(2, 1) = 0.002;
  matrix_h(3, 1) = 0.00138032073697570;
  matrix_h(4, 1) = 0.002;
  matrix_h(5, 1) = 0.00138032073697570;
  matrix_h(6, 1) = 0.000861397758772237;
  matrix_h(7, 1) = 0.000965167479061995;
  matrix_h(8, 1) = 0.000861397758772237;

  matrix_h(0, 2) = 0.000965167479061995;
  matrix_h(1, 2) = 0.00200000000000000;
  matrix_h(2, 2) = 0.00381971863420549;
  matrix_h(3, 2) = 0.000861397758772237;
  matrix_h(4, 2) = 0.00138032073697570;
  matrix_h(5, 2) = 0.002;
  matrix_h(6, 2) = 0.000678804493543927;
  matrix_h(7, 2) = 0.000861397758772237;
  matrix_h(8, 2) = 0.000965167479061995;

  matrix_h(0, 3) = 0.002;
  matrix_h(1, 3) = 0.00138032073697570;
  matrix_h(2, 3) = 0.000861397758772237;
  matrix_h(3, 3) = 0.00381971863420549;
  matrix_h(4, 3) = 0.002;
  matrix_h(5, 3) = 0.000965167479061995;
  matrix_h(6, 3) = 0.002;
  matrix_h(7, 3) = 0.00138032073697570;
  matrix_h(8, 3) = 0.000861397758772237;

  matrix_h(0, 4) = 0.00138032073697570;
  matrix_h(1, 4) = 0.002;
  matrix_h(2, 4) = 0.00138032073697570;
  matrix_h(3, 4) = 0.002;
  matrix_h(4, 4) = 0.00381971863420549;
  matrix_h(5, 4) = 0.002;
  matrix_h(6, 4) = 0.00138032073697570;
  matrix_h(7, 4) = 0.002;
  matrix_h(8, 4) = 0.00138032073697570;

  matrix_h(0, 5) = 0.000861397758772237;
  matrix_h(1, 5) = 0.00138032073697570;
  matrix_h(2, 5) = 0.002;
  matrix_h(3, 5) = 0.000965167479061995;
  matrix_h(4, 5) = 0.002;
  matrix_h(5, 5) = 0.00381971863420549;
  matrix_h(6, 5) = 0.000861397758772237;
  matrix_h(7, 5) = 0.00138032073697570;
  matrix_h(8, 5) = 0.00200000000000000;

  matrix_h(0, 6) = 0.000965167479061995;
  matrix_h(1, 6) = 0.000861397758772237;
  matrix_h(2, 6) = 0.000678804493543927;
  matrix_h(3, 6) = 0.002;
  matrix_h(4, 6) = 0.00138032073697570;
  matrix_h(5, 6) = 0.000861397758772237;
  matrix_h(6, 6) = 0.00381971863420549;
  matrix_h(7, 6) = 0.002;
  matrix_h(8, 6) = 0.000965167479061995;

  matrix_h(0, 7) = 0.000861397758772237;
  matrix_h(1, 7) = 0.000965167479061995;
  matrix_h(2, 7) = 0.000861397758772237;
  matrix_h(3, 7) = 0.00138032073697570;
  matrix_h(4, 7) = 0.002;
  matrix_h(5, 7) = 0.00138032073697570;
  matrix_h(6, 7) = 0.002;
  matrix_h(7, 7) = 0.00381971863420549;
  matrix_h(8, 7) = 0.002;

  matrix_h(0, 8) = 0.000678804493543927;
  matrix_h(1, 8) = 0.000861397758772237;
  matrix_h(2, 8) = 0.000965167479061995;
  matrix_h(3, 8) = 0.000861397758772237;
  matrix_h(4, 8) = 0.00138032073697570;
  matrix_h(5, 8) = 0.002;
  matrix_h(6, 8) = 0.000965167479061995;
  matrix_h(7, 8) = 0.00200000000000000;
  matrix_h(8, 8) = 0.00381971863420549;

  MIRCO::ViewVector_h b0_h("b0_h", 9);
  b0_h(0) = 1357.42803841637;
  b0_h(1) = 1330.45347724905;
  b0_h(2) = 1357.42803841637;
  b0_h(3) = 1353.38917789346;
  b0_h(4) = 1376.01952100849;
  b0_h(5) = 1350.64903939096;
  b0_h(6) = 1357.42803841637;
  b0_h(7) = 1411.27792115093;
  b0_h(8) = 1357.42803841637;

  MIRCO::ViewVector_h p0_h("p0_h", 9);
  p0_h(0) = 161643.031767534;
  p0_h(1) = 43277.7916790043;
  p0_h(2) = 162132.580424512;
  p0_h(3) = 54559.4126169668;
  p0_h(4) = 10518.1563067329;
  p0_h(5) = 53208.9583111822;
  p0_h(6) = 147203.068996657;
  p0_h(7) = 83111.4417592719;
  p0_h(8) = 147692.617653634;

  MIRCO::ViewMatrix_d matrix_d =
      Kokkos::create_mirror_view_and_copy(MIRCO::ExecSpace_Default_t(), matrix_h);
  MIRCO::ViewVector_d p0_d =
      Kokkos::create_mirror_view_and_copy(MIRCO::ExecSpace_Default_t(), p0_h);
  MIRCO::ViewVector_d b0_d =
      Kokkos::create_mirror_view_and_copy(MIRCO::ExecSpace_Default_t(), b0_h);

  MIRCO::ViewVectorInt_d activeSet0_d("activeSet0_d", 9);
  Kokkos::parallel_for(9, KOKKOS_LAMBDA(const int i) { activeSet0_d(i) = i; });

  MIRCO::ViewVector_d pf_d;
  MIRCO::ViewVectorInt_d activeSetf_d;

  MIRCO::nonlinearSolve(pf_d, activeSetf_d, p0_d, activeSet0_d, matrix_d, b0_d);

  Kokkos::deep_copy(p0_h, p0_d);

  EXPECT_NEAR(p0_h(0), 163213.374921086, 1e-06);
  EXPECT_NEAR(p0_h(1), 43877.9231473546, 1e-06);
  EXPECT_NEAR(p0_h(2), 163702.923578063, 1e-06);
  EXPECT_NEAR(p0_h(3), 55159.5440853170, 1e-06);
  EXPECT_NEAR(p0_h(4), 10542.1713862417, 1e-06);
  EXPECT_NEAR(p0_h(5), 53809.0897795325, 1e-06);
  EXPECT_NEAR(p0_h(6), 148773.412150208, 1e-06);
  EXPECT_NEAR(p0_h(7), 83711.5732276221, 1e-06);
  EXPECT_NEAR(p0_h(8), 149262.960807186, 1e-06);
}

TEST(FilesystemUtils, createrelativepath)
{
  std::string targetfilename = "input.dat";
  std::string sourcefilename = "../inputfiles/sourceinput.json";
  MIRCO::Utils::changeRelativePath(targetfilename, sourcefilename);
  EXPECT_EQ(targetfilename, "../inputfiles/input.dat");
}

TEST(FilesystemUtils, keepabsolutepath)
{
  std::string targetfilename = "/root_dir/home/user/Input/input.dat";
  std::string sourcefilename = "../inputfiles/sourceinput.json";
  MIRCO::Utils::changeRelativePath(targetfilename, sourcefilename);
  EXPECT_EQ(targetfilename, "/root_dir/home/user/Input/input.dat");
}

TEST(topology, RMG)
{
  int Resolution = 2;
  float HurstExponent = 0.1;
  bool RandomSeedFlag = false;
  int RandomGeneratorSeed = 95;
  double InitialTopologyStdDeviation = 20.0;

  MIRCO::ViewMatrix_h outsurf_h = MIRCO::CreateRmgSurface(
      Resolution, InitialTopologyStdDeviation, HurstExponent, RandomSeedFlag, RandomGeneratorSeed);

  EXPECT_NEAR(outsurf_h(0, 0), 23.5435469989256, 1e-06);
  EXPECT_NEAR(outsurf_h(0, 1), 30.2624522170979, 1e-06);
  EXPECT_NEAR(outsurf_h(0, 2), 69.5813622417479, 1e-06);
  EXPECT_NEAR(outsurf_h(0, 3), 43.5026425381265, 1e-06);
  EXPECT_NEAR(outsurf_h(0, 4), 23.5435469989256, 1e-06);
  EXPECT_NEAR(outsurf_h(1, 0), 68.8507553267314, 1e-06);
  EXPECT_NEAR(outsurf_h(1, 1), 73.8350740079714, 1e-06);
  EXPECT_NEAR(outsurf_h(1, 2), 77.9927972851754, 1e-06);
  EXPECT_NEAR(outsurf_h(1, 3), 35.2927793006724, 1e-06);
  EXPECT_NEAR(outsurf_h(1, 4), 22.6620325442329, 1e-06);
  EXPECT_NEAR(outsurf_h(2, 0), 39.1583562054882, 1e-06);
  EXPECT_NEAR(outsurf_h(2, 1), 19.2247183888878, 1e-06);
  EXPECT_NEAR(outsurf_h(2, 2), 79.1711886771701, 1e-06);
  EXPECT_NEAR(outsurf_h(2, 3), 5.66729306836534, 1e-06);
  EXPECT_NEAR(outsurf_h(2, 4), 41.3691438722521, 1e-06);
  EXPECT_NEAR(outsurf_h(3, 0), 59.1811726494348, 1e-06);
  EXPECT_NEAR(outsurf_h(3, 1), 21.2400598989696, 1e-06);
  EXPECT_NEAR(outsurf_h(3, 2), 54.6656122080671, 1e-06);
  EXPECT_NEAR(outsurf_h(3, 3), 28.0246974768169, 1e-06);
  EXPECT_NEAR(outsurf_h(3, 4), 6.72730409669533, 1e-06);
  EXPECT_NEAR(outsurf_h(4, 0), 23.5435469989256, 1e-06);
  EXPECT_NEAR(outsurf_h(4, 1), 0, 1e-03);
  EXPECT_NEAR(outsurf_h(4, 2), 30.6777944575233, 1e-06);
  EXPECT_NEAR(outsurf_h(4, 3), 35.2191824993355, 1e-06);
  EXPECT_NEAR(outsurf_h(4, 4), 23.5435469989256, 1e-06);
}
TEST(topology, readFromFile)
{
  int N;
  std::string topologyFilePath = "test/data/topologyN5.dat";
  MIRCO::ViewMatrix_h outsurf_h = MIRCO::CreateSurfaceFromFile(topologyFilePath, N);

  EXPECT_EQ(outsurf_h.extent(0), 5);
  EXPECT_EQ(outsurf_h.extent(1), 5);
  EXPECT_NEAR(outsurf_h(0, 0), 5.7299175e+01, 1e-06);
  EXPECT_NEAR(outsurf_h(4, 3), 9.8243100e+01, 1e-06);
}

TEST(inputParameters, yaml_rmg)
{
  std::string inputFilePath = "test/data/input_res2_rmg.yaml";
  MIRCO::InputParameters inputParams(inputFilePath);

  MIRCO::ViewMatrix_d topology_d = inputParams.topology_d;

  EXPECT_EQ(topology_d.extent(0), 5);
  EXPECT_EQ(topology_d.extent(1), 5);

  EXPECT_EQ(inputParams.max_iteration, 100);
  EXPECT_NEAR(inputParams.tolerance, 0.01, 1e-06);

  EXPECT_NEAR(inputParams.grid_size, 200, 1e-04);
  EXPECT_NEAR(inputParams.composite_youngs, 0.549451, 1e-04);
}
TEST(inputParameters, yaml_dat)
{
  std::string inputFilePath = "test/data/input_withDat.yaml";
  MIRCO::InputParameters inputParams(inputFilePath);

  MIRCO::ViewMatrix_d topology_d = inputParams.topology_d;

  EXPECT_EQ(topology_d.extent(0), 5);
  EXPECT_EQ(topology_d.extent(1), 5);

  MIRCO::ViewMatrix_h topology_h =
      Kokkos::create_mirror_view_and_copy(MIRCO::ExecSpace_DefaultHost_t(), topology_d);
  EXPECT_NEAR(topology_h(0, 0), 5.7299175e+01, 1e-06);
  EXPECT_NEAR(topology_h(4, 3), 9.8243100e+01, 1e-06);

  EXPECT_EQ(inputParams.max_iteration, 100);
  EXPECT_NEAR(inputParams.tolerance, 0.01, 1e-06);

  EXPECT_NEAR(inputParams.grid_size, 200, 1e-04);
  EXPECT_NEAR(inputParams.composite_youngs, 0.549451, 1e-04);
}
TEST(inputParameters, directInput_rmg)
{
  MIRCO::InputParameters inputParams(
      1.0, 1.0, 0.2, 0.2, 0.005, 10.0, 1000, 2, 15.0, 0.15, false, 46, 100, false, true);

  MIRCO::ViewMatrix_d topology_d = inputParams.topology_d;

  EXPECT_EQ(topology_d.extent(0), 5);
  EXPECT_EQ(topology_d.extent(1), 5);

  EXPECT_EQ(inputParams.max_iteration, 100);
  EXPECT_NEAR(inputParams.tolerance, 0.005, 1e-06);

  EXPECT_NEAR(inputParams.grid_size, 200, 1e-04);
  EXPECT_EQ(inputParams.N, 5);
}
TEST(inputParameters, directInput_dat)
{
  std::string topologyFilePath = "test/data/topologyN5.dat";
  MIRCO::InputParameters inputParams(
      1.0, 1.0, 0.2, 0.2, 0.005, 10.0, 1000, topologyFilePath, 100, false, false);

  MIRCO::ViewMatrix_d topology_d = inputParams.topology_d;

  EXPECT_EQ(topology_d.extent(0), 5);
  EXPECT_EQ(topology_d.extent(1), 5);

  MIRCO::ViewMatrix_h topology_h =
      Kokkos::create_mirror_view_and_copy(MIRCO::ExecSpace_DefaultHost_t(), topology_d);
  EXPECT_NEAR(topology_h(0, 0), 5.7299175e+01, 1e-06);
  EXPECT_NEAR(topology_h(4, 3), 9.8243100e+01, 1e-06);

  EXPECT_EQ(inputParams.max_iteration, 100);
  EXPECT_NEAR(inputParams.tolerance, 0.005, 1e-06);

  EXPECT_NEAR(inputParams.grid_size, 200, 1e-04);
  EXPECT_EQ(inputParams.N, 5);
}

TEST(warmstarting, warmstart)
{
  using ViewVectorInt_h = Kokkos::View<int*, Kokkos::LayoutLeft, MIRCO::ExecSpace_DefaultHost_t>;
  ViewVectorInt_h activeSet0_h("activeSet0_h", 3);
  ViewVectorInt_h activeSetf_h("activeSetf_h", 2);
  MIRCO::ViewVector_h pf_h("", 2);

  activeSet0_h(0) = 12;
  activeSet0_h(1) = 34;
  activeSet0_h(2) = 56;

  activeSetf_h(0) = 12;
  activeSetf_h(1) = 56;

  pf_h(0) = 10;
  pf_h(1) = 30;

  MIRCO::ViewVectorInt_d activeSet0_d =
      Kokkos::create_mirror_view_and_copy(MIRCO::MemorySpace_ofDefaultExec_t(), activeSet0_h);
  MIRCO::ViewVectorInt_d activeSetf_d =
      Kokkos::create_mirror_view_and_copy(MIRCO::MemorySpace_ofDefaultExec_t(), activeSetf_h);
  MIRCO::ViewVector_d pf_d =
      Kokkos::create_mirror_view_and_copy(MIRCO::MemorySpace_ofDefaultExec_t(), pf_h);

  MIRCO::ViewVector_d p0_d = MIRCO::Warmstart(activeSet0_d, activeSetf_d, pf_d);
  MIRCO::ViewVector_d p0_h = Kokkos::create_mirror_view_and_copy(MIRCO::MemorySpace_Host_t(), p0_d);

  EXPECT_EQ(p0_h(0), 10);
  EXPECT_EQ(p0_h(1), 0);
  EXPECT_EQ(p0_h(2), 30);
}

TEST(warmstarting, warmstart2)
{
  using ViewVectorInt_h = Kokkos::View<int*, Kokkos::LayoutLeft, MIRCO::ExecSpace_DefaultHost_t>;
  ViewVectorInt_h activeSet0_h("activeSet0_h", 3);
  ViewVectorInt_h activeSetf_h("activeSetf_h", 4);
  MIRCO::ViewVector_h pf_h("", 4);

  activeSet0_h(0) = 12;
  activeSet0_h(1) = 34;
  activeSet0_h(2) = 56;

  activeSetf_h(0) = 12;
  activeSetf_h(1) = 78;
  activeSetf_h(2) = 910;
  activeSetf_h(3) = 56;

  pf_h(0) = 10;
  pf_h(1) = 50;
  pf_h(2) = 70;
  pf_h(3) = 30;

  MIRCO::ViewVectorInt_d activeSet0_d =
      Kokkos::create_mirror_view_and_copy(MIRCO::MemorySpace_ofDefaultExec_t(), activeSet0_h);
  MIRCO::ViewVectorInt_d activeSetf_d =
      Kokkos::create_mirror_view_and_copy(MIRCO::MemorySpace_ofDefaultExec_t(), activeSetf_h);
  MIRCO::ViewVector_d pf_d =
      Kokkos::create_mirror_view_and_copy(MIRCO::MemorySpace_ofDefaultExec_t(), pf_h);

  MIRCO::ViewVector_d p0_d = MIRCO::Warmstart(activeSet0_d, activeSetf_d, pf_d);
  MIRCO::ViewVector_d p0_h = Kokkos::create_mirror_view_and_copy(MIRCO::MemorySpace_Host_t(), p0_d);

  EXPECT_EQ(p0_h(0), 10);
  EXPECT_EQ(p0_h(1), 0);
  EXPECT_EQ(p0_h(2), 30);
}

int main(int argc, char** argv)
{
  Kokkos::initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
