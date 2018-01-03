
AC_INIT([sss],[0.1],[douglas.baxter@pnl.gov])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CC
AC_PROG_CPP

case "$CC" in
	mpi*)
	SERIAL_EXECS=
	SPEC_LIB_EXEC=
	MPI_EXECS=sss_p
	TIMING_CONFIG="no_"
	MMAP_CONFIG="no_"
	CFLAGS="$CFLAGS -DLINUX -DMPI"
	;;
	*)
	SERIAL_EXECS="boltzmann boltzmann_boot_test rxn_map libboltzmann.a"
	MPI_EXECS=
	TIMING_CONFIG="no_"
	MMAP_CONFIG="no_"
	CFLAGS="-O2 -g -DLINUX"
	;;
esac

AC_SUBST(MPI_EXECS)
AC_SUBST(SERIAL_EXECS)
AC_SUBST(TIMING_CONFIG)
AC_SUBST(MMAP_CONFIG)
AC_OUTPUT([Makefile])
