#ifndef SRC_KOKKOSTYPES_KOKKOS_H_
#define SRC_KOKKOSTYPES_KOKKOS_H_

#include <Kokkos_Core.hpp>

// This file defines some commonly used Kokkos type aliases. More aliases, as well as Kokkos-related
// macros, utilities, etc., can be added as needed
namespace MIRCO
{
  using ExecSpace_DefaultHost_t = Kokkos::DefaultHostExecutionSpace;
  using ExecSpace_Default_t = Kokkos::DefaultExecutionSpace;

  using MemorySpace_Host_t = Kokkos::HostSpace;
  using MemorySpace_ofDefaultExec_t = ExecSpace_Default_t::memory_space;

  using Device_Host_t = Kokkos::Device<ExecSpace_DefaultHost_t, MemorySpace_Host_t>;
  using Device_Default_t = Kokkos::Device<ExecSpace_Default_t, MemorySpace_ofDefaultExec_t>;

  // Note: LayoutLeft is column major and is necessary for KokkosLapack::gesv() in
  // MIRCO::nonlinearSolve(). We have to specify it for views of any rank, even if it is
  // physically meaningless for rank < 2.
  using ViewVector_h = Kokkos::View<double*, Kokkos::LayoutLeft, Device_Host_t>;
  using ViewMatrix_h = Kokkos::View<double**, Kokkos::LayoutLeft, Device_Host_t>;

  using ViewVector_d = Kokkos::View<double*, Kokkos::LayoutLeft, Device_Default_t>;
  using ViewMatrix_d = Kokkos::View<double**, Kokkos::LayoutLeft, Device_Default_t>;

  using ViewScalarInt_d = Kokkos::View<int, Kokkos::LayoutLeft, Device_Default_t>;
  using ViewVectorInt_d = Kokkos::View<int*, Kokkos::LayoutLeft, Device_Default_t>;
}  // namespace MIRCO

#endif  // SRC_KOKKOSTYPES_KOKKOS_H_
