#!/bin/bash -e

for ii in $tool_list ; do
    load "${BUILDDIR}/setup_${ii}"
done

# ------------------------------------------------------------------------
# Install or compile packages using newly installed tools
# ------------------------------------------------------------------------

# setup compiler flags, leading to nice stack traces on crashes but
# still optimised
CFLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native $TSANFLAGS"
FFLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native $TSANFLAGS"
F77FLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native $TSANFLAGS"
F90FLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native $TSANFLAGS"
FCFLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native $TSANFLAGS"
CXXFLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native $TSANFLAGS"

export CFLAGS=$(allowed_gcc_flags $CFLAGS)
export FFLAGS=$(allowed_gfortran_flags $FFLAGS)
export F77FLAGS=$(allowed_gfortran_flags $F77FLAGS)
export F90FLAGS=$(allowed_gfortran_flags $F90FLAGS)
export FCFLAGS=$(allowed_gfortran_flags $FCFLAGS)
export CXXFLAGS=$(allowed_gxx_flags $CXXFLAGS)

export LDFLAGS="$TSANFLAGS"

# get system arch information using OpenBLAS prebuild
"${SCRIPTDIR}"/get_openblas_arch.sh; load "${BUILDDIR}/openblas_arch"

# MPI libraries
time_start=`date +%s`
case "$MPI_MODE" in
    mpich)
        "${SCRIPTDIR}"/install_mpich.sh "${with_mpich}"; load "${BUILDDIR}/setup_mpich"
        ;;
    openmpi)
        "${SCRIPTDIR}"/install_openmpi.sh "${with_openmpi}"; load "${BUILDDIR}/setup_openmpi"
        ;;
esac
time_stop=`date +%s`
printf "Step took %0.2f seconds.\n" $((time_stop-time_start))

#EOF