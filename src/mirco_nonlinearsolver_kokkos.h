#ifndef SRC_NONLINEARSOLVER_KOKKOS_H_
#define SRC_NONLINEARSOLVER_KOKKOS_H_

#include "mirco_kokkostypes_kokkos.h"

namespace MIRCO
{
  /**
   * @brief Solve the non-linear set of equations using the Non-Negative Least Squares (NNLS)
   * algorithm (Algorithm 3) from (Bemporad & Paggi, 2015)
   * https://doi.org/10.1016/j.ijsolstr.2015.06.005
   *
   * @param[in] matrix_d Influence coefficient matrix (Discrete version of Green Function)
   * @param[in] b0_d Indentation value of the half space at the predicted points of contact
   * @param[in, out] p_d contact forces vector
   * @param[out] activeSetSize size of the active set
   * @param[in] nnlstol tolerance of the nonlinear solver; \epsilon in (Bemporad & Paggi, 2015)
   * @param[in] maxiter maximum number of total iterations of the innermost loop of the nonlinear
   * solver
   */
  void nonlinearSolve(const ViewMatrix_d matrix_d, const ViewVector_d b0_d, ViewVector_d& p_d,
      int& activeSetSize, double nnlstol = 1.0e-08, int maxiter = 10000);
}  // namespace MIRCO

#endif  // SRC_NONLINEARSOLVER_KOKKOS_H_
