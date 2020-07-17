
/* #undef HAVE_THYRA_EPETRA */

/* #undef HAVE_THYRA_EPETRAEXT */

/* #undef HAVE_THYRA_TPETRA */

/* #undef HAVE_THYRA_TEUCHOS_BLASFLOAT */

#define HAVE_THYRA_FLOAT

#define HAVE_THYRA_COMPLEX

/* #undef HAVE_THYRA_DEBUG */

#define HAVE_THYRA_EXPLICIT_INSTANTIATION

#define HAVE_THYRA_ME_POLYNOMIAL

#ifndef THYRA_FUNC_TIME_MONITOR
#  define THYRA_TEUCHOS_TIME_MONITOR
#  define THYRA_FUNC_TIME_MONITOR(FUNCNAME) \
     TEUCHOS_FUNC_TIME_MONITOR_DIFF(FUNCNAME, THYRA)
#  define THYRA_FUNC_TIME_MONITOR_DIFF(FUNCNAME, DIFF) \
     TEUCHOS_FUNC_TIME_MONITOR_DIFF(FUNCNAME, DIFF)
#endif


#ifndef THYRA_DEPRECATED
#  if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define THYRA_DEPRECATED  __attribute__((__deprecated__))
#  else
#    define THYRA_DEPRECATED
#  endif
#endif

#ifndef THYRA_DEPRECATED_MSG
#  if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 5))
#    define THYRA_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__ (#MSG) ))
#  elif (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define THYRA_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__))
#  else
#    define THYRA_DEPRECATED_MSG(MSG)
#  endif
#endif

