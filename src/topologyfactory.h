#ifndef SRC_TOPOLOGYFACTORY_H_
#define SRC_TOPOLOGYFACTORY_H_

#include <memory>
#include <string>

// forward declaration
class TopologyGeneration;

void CreateSurfaceObject(int n, double Hurst, bool rand_seed_flag, std::string zfilePath,
    bool rmg_flag, std::shared_ptr<TopologyGeneration>& surfacegenerator);

#endif  // SRC_TOPOLOGYFACTORY_H_
