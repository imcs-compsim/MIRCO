#include "topologyfactory.h"
#include <memory>
#include <string>
#include "topology.h"

void CreateSurfaceObject(int n, double Hurst, bool rand_seed_flag, std::string zfilePath,
    bool rmg_flag, std::shared_ptr<TopologyGeneration>& surfacegenerator)
{
  if (rmg_flag)
  {
    surfacegenerator = std::shared_ptr<Rmg>(new Rmg(n, Hurst, rand_seed_flag));
  }
  else
  {
    surfacegenerator = std::shared_ptr<ReadFile>(new ReadFile(n, zfilePath));
  }
}
