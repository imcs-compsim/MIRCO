#ifndef SRC_WARMSTART_KOKKOS_H_
#define SRC_WARMSTART_KOKKOS_H_

#include <Kokkos_Core.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
// tmp
#include "tmpHelpers/Timer.hpp"
#include "tmpHelpers/kokkosIntegration.hpp"


namespace MIRCO
{
  /**
   * @brief This function is used to determine the nodes which were in contact in the last iteration
   * from the current contact set. This helps in making an initial guess of the nodes in contact in
   * the current iteration and speeds up the computation.
   *
   * @param xv0[in] x-coordinates of the points in contact in the previous iteration.
   * @param yv0[in] y-coordinates of the points in contact in the previous iteration.
   * @param xvf[in] x-coordinates of the points in contact in the previous iteration.
   * @param yvf[in] y-coordinates of the points in contact in the previous iteration.
   * @param pf[in] Contact force at (xvf,yvf) predicted in the previous iteration.
   *
   * @return p0 vector of contact forces predicted in the previous iteration but are a part of
   * currect predicted contact set.
   */
  ViewVector_h Warmstart(const ViewVector_h& xv0, const ViewVector_h& yv0, const ViewVector_h& xvf,
      const ViewVector_h& yvf, const ViewVector_h& pf);

  /**
   * @brief Overload to use on default device.
   */
  // ViewVector_d Warmstart(const ViewVector_d xv0, const ViewVector_d yv0, const ViewVector_d xvf,
  //   const ViewVector_d yvf, const ViewVector_d pf);
}  // namespace MIRCO

#endif  // SRC_WARMSTART_KOKKOS_H_
