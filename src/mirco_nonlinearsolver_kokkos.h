#ifndef SRC_NONLINEARSOLVER_H_
#define SRC_NONLINEARSOLVER_H_

#include <Kokkos_Core.hpp>
// tmp
#include "tmpHelpers/Timer.hpp"
#include "tmpHelpers/kokkosIntegration.hpp"

namespace MIRCO
{
  namespace NonLinearSolver
  {
    /**
     * @brief Solve the non-linear set of equations using Non-Negative Least Squares (NNLS) method.
     *
     * @param[in] matrix Influence coefficient matrix (Discrete version of Green Function)
     * @param[in] b0 Indentation value of the half space at the predicted points of contact.
     * @param[in] y0 contact forces vector _{l(i,j)} predicted in the previous iteration but are a
     * part of currect predicted contact set.
     * @param[in] w Gap between the point on the topology and the half space
     *
     * @return y Solution containing force
     */

    // # _h for now; we will see whether we can just keep on device more
    ViewVector_h Solve(const ViewMatrix_h& matrix, const ViewVector_h& b0, const ViewVector_h& y0,
        ViewMatrix_h& w);
  };  // namespace NonLinearSolver
}  // namespace MIRCO

#endif  // SRC_NONLINEARSOLVER_H_
