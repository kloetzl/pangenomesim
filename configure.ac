AC_INIT([pangenomesim], [0.2])
AM_INIT_AUTOMAKE([-Wall foreign])

AC_PROG_CXX
AC_LANG(C++)

AC_CONFIG_MACRO_DIR([m4])
AX_CXX_COMPILE_STDCXX(14, [noext], [mandatory])

AC_TYPE_SIZE_T
AC_CHECK_FUNCS([pow])
AC_CHECK_FUNCS([mkdir])
AC_CHECK_HEADER_STDBOOL

AC_CONFIG_HEADERS([src/config.h:src/config.hin])

AC_CONFIG_FILES([
 Makefile
 pangenomesim.1
])
AC_OUTPUT
