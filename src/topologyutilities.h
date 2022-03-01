#ifndef SRC_TOPOLOGYUTILITIES_H_
#define SRC_TOPOLOGYUTILITIES_H_

#include <Epetra_SerialSymDenseMatrix.h>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

void CreateMeshgrid(std::vector<double>& x, int iter, double delta);

// forward declaration
class TopologyGeneration;

void CreateSurfaceObject(int n, double Hurst, bool rand_seed_flag, std::string zfilePath,
    bool rmg_flag, int rmg_seed, std::shared_ptr<TopologyGeneration>& surfacegenerator);

void ComputeMaxAndMean(Epetra_SerialDenseMatrix topology, double& zmax, double& zmean);

#endif  // SRC_TOPOLOGYUTILITIES_H_