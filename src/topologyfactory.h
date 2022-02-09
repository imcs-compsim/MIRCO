#include <memory>
#include <string>

void CreateSurfaceObject(int n, double Hurst, bool rand_seed_flag,
                         std::string zfilePath, bool rmg_flag,
                         std::shared_ptr<TopologyGeneration> &surfacegenerator);