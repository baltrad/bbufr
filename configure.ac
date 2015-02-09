#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(bufr-opera-mf, 1.21)
AC_CONFIG_SRCDIR([lib/bufr.h])
AM_INIT_AUTOMAKE([foreign no-dependencies no-define])
AC_CONFIG_HEADERS([lib/config.h])
AC_CONFIG_FILES([Makefile lib/Makefile tables/Makefile tests/Makefile view/Makefile])

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_RANLIB

# Check for machname, taken from Python configure scripts
# in order to get correct shared object compilation flag
#
AC_SUBST(MACHDEP)
AC_MSG_CHECKING(MACHDEP)
if test -z "$MACHDEP"
then
	ac_sys_system=`uname -s`
	if test "$ac_sys_system" = "AIX" \
	-o "$ac_sys_system" = "UnixWare" -o "$ac_sys_system" = "OpenUNIX"; then
		ac_sys_release=`uname -v`
	else
		ac_sys_release=`uname -r`
	fi
	ac_md_system=`echo $ac_sys_system |
			   tr -d '[/ ]' | tr '[[A-Z]]' '[[a-z]]'`
	ac_md_release=`echo $ac_sys_release |
			   tr -d '[/ ]' | sed 's/^[[A-Z]]\.//' | sed 's/\..*//'`
	MACHDEP="$ac_md_system$ac_md_release"

	case $MACHDEP in
	cygwin*) MACHDEP="cygwin";;
	darwin*) MACHDEP="darwin";;
	atheos*) MACHDEP="atheos";;
        irix646) MACHDEP="irix6";;
	'')	MACHDEP="unknown";;
	esac
fi
AC_MSG_RESULT($MACHDEP)

# Check for compiler options
AC_MSG_CHECKING(for CCSHARED)
if test -z "$CCSHARED"
then
	case $ac_sys_system/$ac_sys_release in
	SunOS*) if test "$GCC" = yes;
		then CCSHARED="-fPIC";
		elif test `uname -p` = sparc;
		then CCSHARED="-xcode=pic32";
		else CCSHARED="-Kpic";
		fi;;
	hp*|HP*) if test "$GCC" = yes;
		 then CCSHARED="-fPIC";
		 else CCSHARED="+z";
		 fi;;
	Linux*|GNU*) CCSHARED="-fPIC";;
	BSD/OS*/4*) CCSHARED="-fpic";;
	FreeBSD*|NetBSD*|OpenBSD*|DragonFly*) CCSHARED="-fPIC";;
	OpenUNIX*|UnixWare*)
		if test "$GCC" = "yes"
		then CCSHARED="-fPIC"
		else CCSHARED="-KPIC"
		fi;;
	SCO_SV*)
		if test "$GCC" = "yes"
		then CCSHARED="-fPIC"
		else CCSHARED="-Kpic -belf"
		fi;;
	IRIX*/6*)  case $CC in
		   *gcc*) CCSHARED="-shared";;
		   *) CCSHARED="";;
		   esac;;
	atheos*) CCSHARED="-fPIC";;
	esac
fi
AC_MSG_RESULT($CCSHARED)

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
#AC_FUNC_MALLOC
AC_FUNC_MEMCMP
#AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset pow strchr])

CFLAGS="$CFLAGS $CCSHARED"
AC_OUTPUT
