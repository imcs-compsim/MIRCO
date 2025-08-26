#ifndef SRC_NONLINEARSOLVER_H_
#define SRC_NONLINEARSOLVER_H_

#include "mirco_kokkostypes.h"

namespace MIRCO
{
  /**
   * @brief Solve the non-linear problem using a Non-Negative Least Squares (NNLS)
   *
   * This implementation follows the Non-Negative Least Squares (NNLS)
   * algorithm (Algorithm 3) from (Bemporad & Paggi, 2015)
   * https://doi.org/10.1016/j.ijsolstr.2015.06.005
   *
   * @param[out] pf final contact forces vector in compact form (only nonzero forces)
   * @param[out] activeSetf final active set at the end of the nonlinear solver
   * @param[in] p full contact forces vector initial guess
   * @param[in] activeSet0 active set initial guess
   * @param[in] matrix Influence coefficient matrix (Discrete version of Green Function)
   * @param[in] b0 Indentation value of the half space at the predicted points of contact
   * @param[in] nnlstol tolerance of the nonlinear solver; \epsilon in (Bemporad & Paggi, 2015)
   * @param[in] maxiter maximum number of total iterations of the innermost loop of the nonlinear
   * solver
   */
  void nonlinearSolve(ViewVector_d& pf, ViewVectorInt_d& activeSetf, ViewVector_d& p,
      const ViewVectorInt_d activeSet0, const ViewMatrix_d matrix, const ViewVector_d b0,
      double nnlstol = 1.0e-08, int maxiter = 10000);
}  // namespace MIRCO

#endif  // SRC_NONLINEARSOLVER_H_
