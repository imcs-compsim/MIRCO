#ifndef RTOP_CONFIG_H
#define RTOP_CONFIG_H

#define HAVE_MPI

/* #undef HAVE_RTOP_DEBUG */

#define HAVE_RTOP_EXPLICIT_INSTANTIATION

#ifndef RTOP_DEPRECATED
#  if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define RTOP_DEPRECATED  __attribute__((__deprecated__))
#  else
#    define RTOP_DEPRECATED
#  endif
#endif

#ifndef RTOP_DEPRECATED_MSG
#  if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 5))
#    define RTOP_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__ (#MSG) ))
#  elif (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define RTOP_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__))
#  else
#    define RTOP_DEPRECATED_MSG(MSG)
#  endif
#endif


#endif
