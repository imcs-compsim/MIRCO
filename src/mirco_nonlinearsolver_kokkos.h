#ifndef SRC_NONLINEARSOLVER_KOKKOS_H_
#define SRC_NONLINEARSOLVER_KOKKOS_H_

namespace MIRCO
{
  #include "mirco_kokkostypes_kokkos.h"
  /**
    * @brief Solve the non-linear set of equations using the Non-Negative Least Squares (NNLS)
    * algorithm (Algorithm 3) from (Bemporad & Paggi, 2015)
    * https://doi.org/10.1016/j.ijsolstr.2015.06.005
    *
    * @param[in] matrix Influence coefficient matrix (Discrete version of Green Function)
    * @param[in] b0 Indentation value of the half space at the predicted points of contact.
    * @param[in, out] p_d contact forces vector _{l(i,j)} predicted in the previous iteration but
    * are a part of currect predicted contact set.
    * @param[out] w Gap between the point on the topology and the half space
    */
  void nonlinearSolve(const ViewMatrix_d matrix, const ViewVector_d b0_d, ViewVector_d& p_d,
      int& activeSetSize, double nnlstol = 1.0e-08, int maxiter = 10000);
}  // namespace MIRCO

#endif  // SRC_NONLINEARSOLVER_KOKKOS_H_
