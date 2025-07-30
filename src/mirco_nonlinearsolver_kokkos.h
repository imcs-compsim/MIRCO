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
     * method from (Lawson and Hanson 1974) https://doi.org/10.1137/1.9781611971217.ch23
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
    ViewVector_d Solve(const ViewMatrix_d& matrix, const ViewVector_d& b0_d, const ViewVector_h& p0, ViewVector_h& w_out, bool returnWRef = true, double nnlstol = 1.0e-08, int maxiter = 1000);
  };  // namespace NonLinearSolver
}  // namespace MIRCO

#endif  // SRC_NONLINEARSOLVER_KOKKOS_H_
