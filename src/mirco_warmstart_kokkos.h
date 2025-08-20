#ifndef SRC_WARMSTART_KOKKOS_H_
#define SRC_WARMSTART_KOKKOS_H_

#include "mirco_kokkostypes_kokkos.h"

namespace MIRCO
{
  /**
   * @brief This function is used to determine the nodes which were in contact in the last iteration
   * from the current contact set. This helps in making an initial guess of the nodes in contact in
   * the current iteration and speeds up the computation.
   *
   * @param[in] activeSet0_d Points predicted to be in contact in the current iteration
   * @param[in] activeSetf_d Points in contact in the previous iteration
   * @param[in] pf_d Contact force vector predicted in the previous iteration
   *
   * @return p0_d vector of contact forces predicted in the previous iteration but which are a part
   * of the currect predicted contact set (warmstart prediction)
   */
  ViewVector_d Warmstart(
      const ViewVector_d& activeSet0_d, const ViewVector_d& activeSetf_d, const ViewVector_d& pf);
}  // namespace MIRCO

#endif  // SRC_WARMSTART_KOKKOS_H_
