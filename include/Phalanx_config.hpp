#ifndef PHALANX_CONFIG_HPP
#define PHALANX_CONFIG_HPP

#define HAVE_MPI

/* #undef PHX_DEBUG */

#define PHX_TEUCHOS_TIME_MONITOR

/* #undef PHX_ETI */

#define PHX_KOKKOS_DEVICE_TYPE_SERIAL
/* #undef PHX_KOKKOS_DEVICE_TYPE_CUDA */
/* #undef PHX_KOKKOS_DEVICE_TYPE_THREAD */
/* #undef PHX_KOKKOS_DEVICE_TYPE_OPENMP */

/* #undef PHX_INDEX_SIZE_TYPE_KOKKOS */
/* #undef PHX_INDEX_SIZE_TYPE_INT */
#define PHX_INDEX_SIZE_TYPE_UINT
/* #undef PHX_INDEX_SIZE_TYPE_LONGINT */
/* #undef PHX_INDEX_SIZE_TYPE_ULONGINT */

/* #undef HAVE_PHALANX_TVMET */

#define Phalanx_ENABLE_Intrepid

#ifndef PHALANX_DEPRECATED
#  if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define PHALANX_DEPRECATED  __attribute__((__deprecated__))
#  else
#    define PHALANX_DEPRECATED
#  endif
#endif

#ifndef PHALANX_DEPRECATED_MSG
#  if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 5))
#    define PHALANX_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__ (#MSG) ))
#  elif (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define PHALANX_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__))
#  else
#    define PHALANX_DEPRECATED_MSG(MSG)
#  endif
#endif


#endif
