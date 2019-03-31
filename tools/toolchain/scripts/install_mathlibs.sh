#!/bin/bash -e

# get system arch information using OpenBLAS prebuild
"${SCRIPTDIR}"/get_openblas_arch.sh; load "${BUILDDIR}/openblas_arch"

# math core libraries, need to use reflapack for valgrind builds, as
# many fast libraries are not necessarily valgrind clean
export REF_MATH_CFLAGS=''
export REF_MATH_LDFLAGS=''
export REF_MATH_LIBS=''
export FAST_MATH_CFLAGS=''
export FAST_MATH_LDFLAGS=''
export FAST_MATH_LIBS=''

time_start=`date +%s`
"${SCRIPTDIR}"/install_reflapack.sh "${with_reflapack}"; load "${BUILDDIR}/setup_reflapack"
time_stop=`date +%s`
printf "Step took %0.2f seconds.\n" $((time_stop-time_start))

time_start=`date +%s`
case "$FAST_MATH_MODE" in
    mkl)
        "${SCRIPTDIR}"/install_mkl.sh "${with_mkl}"; load "${BUILDDIR}/setup_mkl"
        ;;
    acml)
        "${SCRIPTDIR}"/install_acml.sh "${with_acml}"; load "${BUILDDIR}/setup_acml"
        ;;
    openblas)
        "${SCRIPTDIR}"/install_openblas.sh "${with_openblas}"; load "${BUILDDIR}/setup_openblas"
        ;;
    cray)
        # note the space is intentional so that the variable is
        # non-empty and can pass require_env checks
        export FAST_MATH_LDFLAGS="${FAST_MATH_LDFLAGS} "
        export FAST_MATH_LIBS="${FAST_MATH_LIBS} ${CRAY_EXTRA_LIBS}"
        ;;
esac
time_stop=`date +%s`
printf "Step took %0.2f seconds.\n" $((time_stop-time_start))

if [ $ENABLE_VALGRIND = "__TRUE__" ] ; then
    export MATH_CFLAGS="${REF_MATH_CFLAGS}"
    export MATH_LDFLAGS="${REF_MATH_LDFLAGS}"
    export MATH_LIBS="${REF_MATH_LIBS}"
else
    export MATH_CFLAGS="${FAST_MATH_CFLAGS}"
    export MATH_LDFLAGS="${FAST_MATH_LDFLAGS}"
    export MATH_LIBS="${FAST_MATH_LIBS}"
fi

export CP_CFLAGS="${CP_CFLAGS} IF_DEBUG(${REF_MATH_CFLAGS}|IF_VALGRIND(${REF_MATH_CFLAGS}|${FAST_MATH_CFLAGS}))"
export CP_LDFLAGS="${CP_LDFLAGS} IF_DEBUG(${REF_MATH_LDFLAGS}|IF_VALGRIND(${REF_MATH_LDFLAGS}|${FAST_MATH_LDFLAGS}))"
export CP_LIBS="${CP_LIBS} IF_DEBUG(${REF_MATH_LIBS}|IF_VALGRIND(${REF_MATH_LIBS}|${FAST_MATH_LIBS}))"

#EOF