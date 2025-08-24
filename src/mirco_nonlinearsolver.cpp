#include "mirco_nonlinearsolver.h"

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
  void nonlinearSolve(ViewVector_d& pf, ViewVectorInt_d& activeSetf, ViewVector_d& p,
      const ViewVectorInt_d activeSet0, const ViewMatrix_d matrix, const ViewVector_d b0,
      double nnlstol, int maxiter)
  {
    using minloc_t = Kokkos::MinLoc<double, int, MemorySpace_ofDefaultExec_t>;
    using minloc_value_t = typename minloc_t::value_type;
    const std::string kokkosLabelPrefix = "nonlinearSolve(); ";

    constexpr double eps = std::numeric_limits<double>::epsilon();

    const int n0 = b0.extent(0);

    ViewVector_d w(kokkosLabelPrefix + "w", n0);

    ViewVectorInt_d activeInactiveSet(kokkosLabelPrefix + "activeInactiveSet", n0);

    ViewScalarInt_d counterActive(kokkosLabelPrefix + "counterActive");
    Kokkos::deep_copy(counterActive, 0);
    ViewScalarInt_d counterInactive(kokkosLabelPrefix + "counterInactive");
    Kokkos::deep_copy(counterInactive, 0);
    Kokkos::parallel_for(
        n0, KOKKOS_LAMBDA(const int i) {
          if (p(i) >= nnlstol)
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
      Kokkos::parallel_for(n0, KOKKOS_LAMBDA(const int i) { w(i) = -b0(i); });
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
            const double wi = w(i);
            if (wi < ml.val)
            {
              ml.val = wi;
              ml.loc = i;
            }
          },
          minloc_t(minloc_w_i_d));

      minloc_value_t minloc_w_i_h;
      Kokkos::deep_copy(minloc_w_i_h, minloc_w_i_d);

      if (init && (minloc_w_i_h.val >= -nnlstol)) break;

      if (init)
      {
        Kokkos::View<minloc_value_t, MemorySpace_ofDefaultExec_t> minloc_w_iInactive_d(
            kokkosLabelPrefix + "minloc_w_iInactive_d");
        Kokkos::parallel_reduce(
            Kokkos::RangePolicy<ExecSpace_Default_t>(activeSetSize, n0),
            KOKKOS_LAMBDA(const int i, minloc_value_t& ml) {
              const double wi = w(activeInactiveSet(i));
              if (wi < ml.val)
              {
                ml.val = wi;
                ml.loc = i;
              }
            },
            minloc_t(minloc_w_iInactive_d));

        minloc_value_t minloc_w_iInactive_h;
        Kokkos::deep_copy(minloc_w_iInactive_h, minloc_w_iInactive_d);

        swapEntries(activeInactiveSet, minloc_w_iInactive_h.loc, activeSetSize);
        ++activeSetSize;
      }
      else
        init = true;

      while (true)
      {
        ++iter;

        // Compact versions of H and b0, i.e. H_I and \overbar{u}_I in line 6 of Algorithm 3,
        // (Bemporad & Paggi, 2015)
        ViewVector_d b0s_compact(kokkosLabelPrefix + "b0_compact", activeSetSize);
        if (activeSetSize > 1)
        {
          ViewMatrix_d H_compact(kokkosLabelPrefix + "H_compact", activeSetSize, activeSetSize);

          Kokkos::parallel_for(
              activeSetSize, KOKKOS_LAMBDA(const int i) {
                const int row = activeInactiveSet(i);
                b0s_compact(i) = b0(row);
                for (int j = 0; j < activeSetSize; ++j)
                {
                  const int col = activeInactiveSet(j);
                  H_compact(i, j) = matrix(row, col);
                }
              });

          ViewVectorInt_d ipiv(kokkosLabelPrefix + "ipiv", activeSetSize);

          // Solve H_I s_I = b0_I; b0s_compact_d becomes s_I
          KokkosLapack::gesv(H_compact, b0s_compact, ipiv);
        }
        else if (activeSetSize == 1)
        {
          Kokkos::parallel_for(
              1, KOKKOS_LAMBDA(const int) {
                const int ii = activeInactiveSet(0);
                b0s_compact(0) = b0(ii) / matrix(ii, ii);
              });
        }

        bool allGreater = true;
        Kokkos::parallel_reduce(
            activeSetSize,
            KOKKOS_LAMBDA(const int i, bool& rallGreater) {
              const bool lesser = (b0s_compact(i) < -nnlstol);
              rallGreater = rallGreater && !lesser;
            },
            Kokkos::LAnd<bool>(allGreater));
        if (allGreater)
        {
          Kokkos::parallel_for(
              activeSetSize,
              KOKKOS_LAMBDA(const int i) { p(activeInactiveSet(i)) = b0s_compact(i); });

          Kokkos::parallel_for(
              n0, KOKKOS_LAMBDA(const int i) {
                double sum = 0.0;
                for (int j = 0; j < activeSetSize; ++j)
                  sum += matrix(i, activeInactiveSet(j)) * b0s_compact(j);
                w(i) = sum - b0(i);
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
                if (b0s_compact(i) <= 0)
                {
                  const double alphai =
                      p(activeInactiveSet(i)) / (eps + p(activeInactiveSet(i)) - b0s_compact(i));
                  if (alphai < ml.val)
                  {
                    ml.val = alphai;
                    ml.loc = i;
                  }
                }
              },
              minloc_t(minloc_alpha_i_d));

          minloc_value_t minloc_alpha_i_h;
          Kokkos::deep_copy(minloc_alpha_i_h, minloc_alpha_i_d);

          Kokkos::parallel_for(
              activeSetSize, KOKKOS_LAMBDA(const int i) {
                p(activeInactiveSet(i)) =
                    p(activeInactiveSet(i)) +
                    minloc_alpha_i_d().val * (b0s_compact(i) - p(activeInactiveSet(i)));
              });

          if (minloc_alpha_i_h.loc > -1)
          {
            Kokkos::deep_copy(Kokkos::subview(p, activeInactiveSet(minloc_alpha_i_h.loc)), 0.0);
            // Note: `swapEntries()` may not be very efficient on GPU because it does not use
            // `Kokkos::kokkos_swap()`, as that can only be used in parallel lambda
            swapEntries(activeInactiveSet, minloc_alpha_i_h.loc, --activeSetSize);
          }
        }
      }
    }
    // Construct the final active set (the lower half of activeInactiveSet), as well as the compact
    // final pressure vector
    activeSetf = ViewVectorInt_d("activeSetf", activeSetSize);
    pf = ViewVector_d("pf", activeSetSize);
    Kokkos::parallel_for(
        activeSetSize, KOKKOS_LAMBDA(const int i) {
          activeSetf(i) = activeSet0(activeInactiveSet(i));
          pf(i) = p(activeInactiveSet(i));
        });
  }

}  // namespace MIRCO
