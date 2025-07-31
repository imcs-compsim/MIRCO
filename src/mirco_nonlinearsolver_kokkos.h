#ifndef SRC_NONLINEARSOLVER_KOKKOS_H_
#define SRC_NONLINEARSOLVER_KOKKOS_H_

#include <Kokkos_Core.hpp>
// tmp
#include "tmpHelpers/kokkosIntegration.hpp"

namespace MIRCO
{
  namespace NonLinearSolver
  {
    /**
     * @brief Solve the non-linear set of equations using the Non-Negative Least Squares (NNLS)
     * algorithm (Algorithm 3) from (Bemporad & Paggi, 2015) https://doi.org/10.1016/j.ijsolstr.2015.06.005
     *
     * @param[in] matrix Influence coefficient matrix (Discrete version of Green Function)
     * @param[in] b0 Indentation value of the half space at the predicted points of contact.
     * @param[in, out] p_d contact forces vector _{l(i,j)} predicted in the previous iteration but are a
     * part of currect predicted contact set.
     * @param[out] w Gap between the point on the topology and the half space
     */

    // # _h for now; we will see whether we can just keep on device more
    void solve(const ViewMatrix_d matrix, const ViewVector_d b0_d, ViewVector_d& p_d, int& activeSetSize, double nnlstol = 1.0000e-08, int maxiter = 10000);
    
    /**
     * @brief Compact and sort a ViewVectord_d w.r.t. a secondary ViewVectori_d. For example, p_d is in the order of the original positional index, but some of those are 0 because they are inactive; use activeSet to
     *
     * @param[in] matrix Influence coefficient matrix (Discrete version of Green Function)
     * @param[in] b0 Indentation value of the half space at the predicted points of contact.
     * @param[in] y0 contact forces vector _{l(i,j)} predicted in the previous iteration but are a
     * part of currect predicted contact set.
     * @param[out] w Gap between the point on the topology and the half space
     *
     * @return y Solution containing force
     */

    // # _h for now; we will see whether we can just keep on device more
    // [[nodiscard]]
    // ViewVector_d compactAndSort(ViewVector_d& p_d, ViewVectorInt_d& activeInactiveSet, bool returnWRef = true, double nnlstol = 1.0e-08, int maxiter = 1000);
  };  // namespace NonLinearSolver
}  // namespace MIRCO

#endif  // SRC_NONLINEARSOLVER_KOKKOS_H_
