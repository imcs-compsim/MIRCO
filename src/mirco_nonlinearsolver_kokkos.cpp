#include "mirco_nonlinearsolver_kokkos.h"

#include <KokkosLapack_gesv.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <vector>

// tmp
#include "tmpHelpers/Timer.hpp"

void swap_entries(ViewVector_d v, const int i, const int j) {
  if constexpr (HOSTONLY) {
    std::swap(v(i), v(j));
  }
  else {
    auto vi = Kokkos::subview(v, i);
    auto vj = Kokkos::subview(v, j);
    
    ViewScalarDouble_d tmp("tmp");
    
    Kokkos::deep_copy(tmp, vi);
    Kokkos::deep_copy(vi, vj);
    Kokkos::deep_copy(vj, tmp);
  }
}


ViewVector_d MIRCO::NonLinearSolver::Solve(const ViewMatrix_d& matrix, const ViewVector_d& b0_d, const ViewVector_h& p0, ViewVector_h& w_out, bool returnWRef, double nnlstol, int maxiter)
{
  const std::string thisFctName =
      "MIRCO::NonLinearSolver::Solve";
      
  // 2^{-52}; double machine precision, i.e. smallest number such that 1.0 + eps > 1.0
  constexpr double eps = 1.0 / (1ULL << 52);
  // positive infinity for double
  constexpr double infty = 0x7FF0000000000000;
  const int n0 = b0_d.extent(0);
  
  ViewScalarInt_d counter_init_d("counter_init_d");
  ViewVectorInt_d activeInactiveSet("activeInactiveSet", n0);
  Kokkos::parallel_scan("build initial active set", n0,
    KOKKOS_LAMBDA(const int i, int& counterActive, int& counterInactive, const bool final) {
      if (p0(i) > 0.0)
      {
        if (final) activeInactiveSet(counterActive) = i;
        ++counterActive;
      }
      else {
        if(final) activeInactiveSet(n0 - counterInactive) = i;
        ++counterInactive;
      }
      if (final && i == n0 - 1) counter_d_init = counterActive;
    });
  int activeSetSize = Kokkos::create_mirror_view_and_copy(MemorySpace_Host_t, counter_d_init)(0);

  ViewVector_d w_d("w_d", n0);
  
  bool init = false;
  if (activeSetSize == 0)
  {
    init = true;
    
    // Hp - b0 for p = 0
    Kokkos::parallel_for(
        thisFctName + "__firstParallel", Kokkos::RangePolicy<ExecSpace_Default_t>(0, n0),
        KOKKOS_LAMBDA(const int x) { w_d(x) = -b0_d(x); });
  }
  ViewVector_d p_d("p_d", n0);

  // iter is incremented in while(2)
  int iter = 0;
  while (1)
  {
    double min_w_iActive;
    int argmin_w_iActive;
    Kokkos::parallel_reduce(
        "min_w", w_d.extent(0),
        KOKKOS_LAMBDA(const int i, Kokkos::MinLoc<double, int>::value_type& ml) {
          auto wi = w_d(activeInactiveSet(i))
          if (wi < ml.val)
          {
            ml.val = wi;
            ml.loc = i;
          }
        },
        Kokkos::MinLoc<double, int>(min_w_iActive, argmin_w_iActive));

    if (((min_w_iActive >= -epsilon || activeSetSize == n0) && init) || iter >= maxiter)
      break;

    if (init) {
      double min_w_iInactive;
      int argmin_w_iInactive;
      Kokkos::parallel_reduce("min_w", Kokkos::RangePolicy<ExecSpace_Default_t>(activeSetSize, n0),
        KOKKOS_LAMBDA(const int i, Kokkos::MinLoc<double, int>::value_type& ml) {
          const double wi = w_d(activeInactiveSet(i))
          if (wi < ml.val)
          {
            ml.val = wi;
            ml.loc = i;
          }
        },
        Kokkos::MinLoc<double, int>(min_w_iInactive, argmin_w_iInactive));
          
      swap_entries(activeInactiveSet, argmin_w_iInactive, activeSetSize);
      ++activeSetSize;
    }
    else
      init = true;

    while (2)
    {
      ++iter;
          
      ViewMatrix_d H_compact_d("H_compact", activeSetSize, activeSetSize);
      ViewVector_d b0s_compact_d("b0_compact", activeSetSize);

      Kokkos::parallel_for("define compact views", activeSetSize,
        KOKKOS_LAMBDA(const int i) {
          int row = activeInactiveSet(i);
          b0s_compact_d(i) = b0_d(row);
          for (int j = 0; j < activeSetSize; ++j) {
            int col = activeInactiveSet(j);
            H_compact_d(i, j) = H(row, col);
          }
        });

      // Solve H s = b0
      KokkosLapack::gesv(H_compact_d, b0s_compact_d);

      bool allGreater = true;
      Kokkos::parallel_reduce("allGreater", activeSetSize,
        KOKKOS_LAMBDA(const int i, bool& local_allGreater) {
          if (b0s_compact_d(i) < -nnlstol)
            l_allGreater = false;
        },
        allGreater);

      if (allGreater)
      {
        Kokkos::parallel_for("p = s_compact", activeSetSize, KOKKOS_LAMBDA(const int i) {
            p_d(activeSet(i)) = b0s_compact_d(i);
          });
          
        Kokkos::parallel_for("w = H p - b0", n0,
          KOKKOS_LAMBDA(const int i) {
            double sum = 0.0;
            for(int j=0; j<activeSetSize; ++j)
              sum += matrix(i, activeInactiveSet(j)) * b0s_compact_d(j);
            w_d(i) = sum - b0(i);
          });
        
        break;
      }
      else//# Algo3 line(8)-(11)
      {
        // Searching for minimum value with index position
        double min_alpha_i = infty;
        int argmin_alpha_i = -1;
        
        Kokkos::parallel_reduce(
          "min_alpha", activeSetSize,
          KOKKOS_LAMBDA(const int i, Kokkos::MinLoc<double, int>::value_type& ml) {
            if(s0(activeInactiveSet(i)) <= 0) {
              const double alphai = p_d(activeInactiveSet(i)) / (eps + p_d(activeInactiveSet(i)) - s0(activeInactiveSet(i)));
              if (alphai < ml.val)
              {
                ml.val = alphai;
                ml.loc = i;
              }
            }
          },
          Kokkos::MinLoc<double, int>(min_alpha_i, argmin_alpha_i));
          
        Kokkos::parallel_for(
          "assemble reduced system", activeSetSize, KOKKOS_LAMBDA(const int i) {
            p_d(activeInactiveSet(i)) = p_d(activeInactiveSet(i)) + alpha * (s0(activeInactiveSet(i)) - p_d(activeInactiveSet(i)))
          });
          
        if (j > -1)
        {
          swap_entries(activeInactiveSet, j, --activeSetSize);
          Kokkos::deep_copy(Kokkos::subview(p_d, j), 0);
        }
      }
    }
  }
  
  //# TODO: sort the active set region of activeInactiveSet and at the same time sort p_d and w_d (if needed?) in the same way
  //# wait but if all we need is the final total pressure and the contact area, p_d does not need to be ordered, so that it would be most efficient to integrate the last two functions into here. We can also make a switch or overload for whether we need the (ordered)
  
  // p_d = p;
  // u_star = w + \overbar{u}

  return p_d;
}
