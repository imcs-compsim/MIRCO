#ifndef SRC_NONLINEARSOLVER_KOKKOS_H_
#define SRC_NONLINEARSOLVER_KOKKOS_H_

#include "mirco_kokkostypes_kokkos.h"

namespace MIRCO
{
  /**
   * @brief Solve the non-linear problem using a Non-Negative Least Squares (NNLS)
   *
   * This implementation follows the Non-Negative Least Squares (NNLS)
   * algorithm (Algorithm 3) from (Bemporad & Paggi, 2015)
   * https://doi.org/10.1016/j.ijsolstr.2015.06.005
   *
   * @param[out] pf_d final contact forces vector in compact form (only nonzero forces)
   * @param[out] activeSetf_d final active set at the end of the nonlinear solver
   * @param[in] p_d full contact forces vector initial guess
   * @param[in] activeSet0_d active set initial guess
   * @param[in] matrix_d Influence coefficient matrix (Discrete version of Green Function)
   * @param[in] b0_d Indentation value of the half space at the predicted points of contact
   * @param[in] nnlstol tolerance of the nonlinear solver; \epsilon in (Bemporad & Paggi, 2015)
   * @param[in] maxiter maximum number of total iterations of the innermost loop of the nonlinear
   * solver
   */
  void nonlinearSolve(ViewVector_d& pf_d, ViewVectorInt_d& activeSetf_d, ViewVector_d& p_d,
      const ViewVectorInt_d activeSet0_d, const ViewMatrix_d matrix_d, const ViewVector_d b0_d,
      double nnlstol = 1.0e-08, int maxiter = 10000);
}  // namespace MIRCO

#endif  // SRC_NONLINEARSOLVER_KOKKOS_H_
