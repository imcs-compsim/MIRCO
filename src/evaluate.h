#ifndef SRC_EVALUATE_H_
#define SRC_EVALUATE_H_

#include <string>

namespace MIRCO
{
  void Evaluate(double &force, double Delta, double lato, double delta, int resolution,
      double Hurst, bool rand_seed_flag, std::string zfilePath, bool rmg_flag, bool rmg_seed,
      double errf, double to1, int max_iter, double E, bool flagwarm, double k_el);
}


#endif  // SRC_EVALUATE_H_
