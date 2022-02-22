#include "computecontactforceandarea.h"
#include <cmath>
#include <vector>

void ComputeContactForceAndArea(std::vector<double> &force0, std::vector<double> &area0, int &iter,
    double &w_el, double nf, std::vector<double> pf, int k, double delta, double lato, double k_el)
{
  force0.push_back(0);
  double sum = 0;
  iter = ceil(nf);
#pragma omp parallel for schedule(static, 16) reduction(+ : sum)  // Always same workload -> Static!
  for (int i = 0; i < iter; i++)
  {
    sum += pf[i];
  }
  force0[k] += sum;
  area0.push_back(nf * (pow(delta, 2) / pow(lato, 2)) * 100);
  w_el = force0[k] / k_el;
}