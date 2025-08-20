#include "mirco_nonlinearsolver_kokkos.h"

#include <KokkosLapack_gesv.hpp>

namespace
{
  using namespace MIRCO;
  void swapEntries(ViewVectorInt_d v, const int i, const int j)
  {
    if constexpr (std::is_same<MemorySpace_Host_t, MemorySpace_ofDefaultExec_t>())
    {
      std::swap(v(i), v(j));
    }
    else
    {
      auto vi = Kokkos::subview(v, i);
      auto vj = Kokkos::subview(v, j);

      ViewScalarInt_d tmp("swapEntries(); tmp");

      Kokkos::deep_copy(tmp, vi);
      Kokkos::deep_copy(vi, vj);
      Kokkos::deep_copy(vj, tmp);
    }
  }
}  // namespace

namespace MIRCO
{
  void nonlinearSolve(ViewVector_d& p_d, ViewVector_d& activeSet, const ViewVector_d activeSet0_d,
      const ViewMatrix_d matrix_d, const ViewVector_d b0_d, double nnlstol, int maxiter)
  {
    using minloc_t = Kokkos::MinLoc<double, int, MemorySpace_ofDefaultExec_t>;
    using minloc_value_t = typename minloc_t::value_type;
    const std::string kokkosLabelPrefix = "nonlinearSolve(); ";

    // 2^{-52}; double machine precision, i.e. smallest number such that 1.0 + eps > 1.0
    constexpr double eps = std::numeric_limits<double>::epsilon();

    const int n0 = b0_d.extent(0);

    ViewVector_d w_d("w_d", n0);

    ViewVectorInt_d activeInactiveSet(kokkosLabelPrefix + "activeInactiveSet", n0);

    ViewScalarInt_d counterActive(kokkosLabelPrefix + "counterActive");
    Kokkos::deep_copy(counterActive, 0);
    ViewScalarInt_d counterInactive(kokkosLabelPrefix + "counterInactive");
    Kokkos::deep_copy(counterInactive, 0);
    Kokkos::parallel_for(
        n0, KOKKOS_LAMBDA(const int i) {
          if (p_d(i) >= nnlstol)
          {
            // Note: atomic_fetch_add() returns the old value
            activeInactiveSet(Kokkos::atomic_fetch_add(&counterActive(), 1)) = i;
          }
          else
          {
            activeInactiveSet(n0 - 1 - (Kokkos::atomic_fetch_add(&counterInactive(), 1))) = i;
          }
        });
    int activeSetSize;
    Kokkos::deep_copy(activeSetSize, counterActive);

    bool init = false;
    if (activeSetSize == 0)
    {
      init = true;
      // Hp - b0 for p = 0
      Kokkos::parallel_for(n0, KOKKOS_LAMBDA(const int i) { w_d(i) = -b0_d(i); });
    }

    int iter = 0;
    while (true)
    {
      if ((init && (activeSetSize == n0)) || iter >= maxiter) break;
      Kokkos::View<minloc_value_t, MemorySpace_ofDefaultExec_t> minloc_w_i_d(
          kokkosLabelPrefix + "minloc_w_i_d");

      Kokkos::parallel_reduce(
          n0,
          KOKKOS_LAMBDA(const int i, minloc_value_t& ml) {
            const double wi = w_d(i);
            if (wi < ml.val)
            {
              ml.val = wi;
              ml.loc = i;
            }
          },
          minloc_t(minloc_w_i_d));

      minloc_value_t minloc_w_i;
      Kokkos::deep_copy(minloc_w_i, minloc_w_i_d);

      if (init && (minloc_w_i.val >= -nnlstol)) break;

      if (init)
      {
        Kokkos::View<minloc_value_t, MemorySpace_ofDefaultExec_t> minloc_w_iInactive_d(
            kokkosLabelPrefix + "minloc_w_iInactive_d");
        Kokkos::parallel_reduce(
            Kokkos::RangePolicy<ExecSpace_Default_t>(activeSetSize, n0),
            KOKKOS_LAMBDA(const int i, minloc_value_t& ml) {
              const double wi = w_d(activeInactiveSet(i));
              if (wi < ml.val)
              {
                ml.val = wi;
                ml.loc = i;
              }
            },
            minloc_t(minloc_w_iInactive_d));

        minloc_value_t minloc_w_iInactive;
        Kokkos::deep_copy(minloc_w_iInactive, minloc_w_iInactive_d);

        swapEntries(activeInactiveSet, minloc_w_iInactive.loc, activeSetSize);
        ++activeSetSize;
      }
      else
        init = true;

      while (true)
      {
        ++iter;

        // Compact versions of H and b0, i.e. H_I and \overbar{u}_I in line 6 of Algorithm 3,
        // (Bemporad & Paggi, 2015)
        ViewVector_d b0s_compact_d(kokkosLabelPrefix + "b0_compact_d", activeSetSize);
        if (activeSetSize > 1)
        {
          ViewMatrix_d H_compact_d(kokkosLabelPrefix + "H_compact_d", activeSetSize, activeSetSize);

          Kokkos::parallel_for(
              activeSetSize, KOKKOS_LAMBDA(const int i) {
                const int row = activeInactiveSet(i);
                b0s_compact_d(i) = b0_d(row);
                for (int j = 0; j < activeSetSize; ++j)
                {
                  const int col = activeInactiveSet(j);
                  H_compact_d(i, j) = matrix_d(row, col);
                }
              });

          ViewVectorInt_d ipiv_d(kokkosLabelPrefix + "ipiv_d", activeSetSize);

          // Solve H_I s_I = b0_I; b0s_compact_d becomes s_I
          KokkosLapack::gesv(H_compact_d, b0s_compact_d, ipiv_d);
        }
        else if (activeSetSize == 1)
        {
          Kokkos::parallel_for(
              1, KOKKOS_LAMBDA(const int) {
                const int ii = activeInactiveSet(0);
                b0s_compact_d(0) = b0_d(ii) / matrix_d(ii, ii);
              });
        }

        bool allGreater = true;
        Kokkos::parallel_reduce(
            activeSetSize,
            KOKKOS_LAMBDA(const int i, bool& rallGreater) {
              const bool lesser = (b0s_compact_d(i) < -nnlstol);
              rallGreater = rallGreater && !lesser;
            },
            Kokkos::LAnd<bool>(allGreater));
        if (allGreater)
        {
          Kokkos::parallel_for(
              activeSetSize,
              KOKKOS_LAMBDA(const int i) { p_d(activeInactiveSet(i)) = b0s_compact_d(i); });

          Kokkos::parallel_for(
              n0, KOKKOS_LAMBDA(const int i) {
                double sum = 0.0;
                for (int j = 0; j < activeSetSize; ++j)
                  sum += matrix_d(i, activeInactiveSet(j)) * b0s_compact_d(j);
                w_d(i) = sum - b0_d(i);
              });

          break;
        }
        else
        {
          // Find min_i and argmin_i of alpha_i := \frac{p_i}{p_i - s_i}
          Kokkos::View<minloc_value_t, MemorySpace_ofDefaultExec_t> minloc_alpha_i_d(
              kokkosLabelPrefix + "minloc_alpha_i_d");
          Kokkos::parallel_reduce(
              activeSetSize,
              KOKKOS_LAMBDA(const int i, minloc_value_t& ml) {
                if (b0s_compact_d(i) <= 0)
                {
                  const double alphai = p_d(activeInactiveSet(i)) /
                                        (eps + p_d(activeInactiveSet(i)) - b0s_compact_d(i));
                  if (alphai < ml.val)
                  {
                    ml.val = alphai;
                    ml.loc = i;
                  }
                }
              },
              minloc_t(minloc_alpha_i_d));

          minloc_value_t minloc_alpha_i;
          Kokkos::deep_copy(minloc_alpha_i, minloc_alpha_i_d);

          Kokkos::parallel_for(
              activeSetSize, KOKKOS_LAMBDA(const int i) {
                p_d(activeInactiveSet(i)) =
                    p_d(activeInactiveSet(i)) +
                    minloc_alpha_i_d().val * (b0s_compact_d(i) - p_d(activeInactiveSet(i)));
              });

          if (minloc_alpha_i.loc > -1)
          {
            Kokkos::deep_copy(Kokkos::subview(p_d, activeInactiveSet(minloc_alpha_i.loc)), 0.0);
            // Note: `swapEntries()` may not be very efficient on GPU because it does not use
            // `Kokkos::kokkos_swap()`, as that can only be used in parallel lambda
            swapEntries(activeInactiveSet, minloc_alpha_i.loc, --activeSetSize);
          }
        }
      }
    }
    // Construct the final active set (the lower half of activeInactiveSet)
    // #
    // # WAIT NO THIS IS WRONG; THIS TAKES INDEX FROM SUBSET OF THE GIVEN ACTIVE SET HERE, BUT WE
    // NEED TO USE
    activeSet = ViewVector_d(kokkosLabelPrefix + "activeSet", activeSetSize);
    Kokkos::parallel_for(
        activeSetSize, KOKKOS_LAMBDA(const int i) {
          activeSet(activeSetSize) = activeSet0_d(activeInactiveSet(activeSetSize));
        });
  }

}  // namespace MIRCO
