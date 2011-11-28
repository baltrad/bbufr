#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(bufr-opera-mf, 1.18)
AC_CONFIG_SRCDIR([lib/bufr.h])
AM_INIT_AUTOMAKE([foreign no-dependencies no-define])
AC_CONFIG_HEADERS([lib/config.h])
AC_CONFIG_FILES([Makefile lib/Makefile tables/Makefile tests/Makefile view/Makefile])

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_RANLIB

# Checks for libraries.
AC_CHECK_LIB([m], [pow])
AC_CHECK_LIB([proj], [pj_init_plus],, 
  [AC_MSG_ERROR([cannot find PROJ library, fatal])])
AC_CHECK_LIB([z], [compress2],, 
  [AC_MSG_ERROR([cannot find zlib library, fatal])])

AC_CHECK_LIB([hdf5], [H5Acreate],, 
  [AC_CHECK_LIB([hdf5], [H5Acreate2],
[AC_DEFINE(H5_USE_16_API,1,[Define to 1 if you have the HDF5 1.8 library.])
AC_DEFINE(HAVE_LIBHDF5,1,[Define to 1 if you have the HDF5 library (-lhdf5).])
LIBS="-lhdf5 $LIBS"],
[AC_MSG_ERROR([cannot find HDF5 library, fatal])])])

AC_CHECK_LIB([png], [png_write_info],, 
  [AC_MSG_WARN([cannot find PNG library, view directory programs won't be built])])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([memory.h stdlib.h string.h proj_api.h zlib.h hdf5.h],,
  [AC_MSG_ERROR([cannot find some system header file, fatal])])
AC_CHECK_HEADERS([png.h],,
  [AC_MSG_WARN([cannot find PNG header file, view directory programs won't be built])])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset pow strchr])

AC_OUTPUT
