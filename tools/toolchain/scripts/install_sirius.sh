#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

sirius_ver=${sirius_ver:-6.1.2}
source "${SCRIPT_DIR}"/step_begin.sh
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh

with_sirius=${1:-__INSTALL__}

if [ "$MPI_MODE" = "no" ] && [ "$ENABLE_OMP" = "__FALSE__" ] ; then
    report_warning $LINENO "MPI and OpenMP are disabled, skipping sirius installation"
    echo 'with_sirius="__FALSE__"' >> ${BUILDDIR}/setup_sirius
    exit 0
fi

[ -f "${BUILDDIR}/setup_sirius" ] && rm "${BUILDDIR}/setup_sirius"

SIRIUS_CFLAGS=''
SIRIUS_LDFLAGS=''
SIRIUS_LIBS=''

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
case "$with_sirius" in
    __DONTUSE__)
    ;;
    __INSTALL__)
        echo "==================== Installing SIRIUS ===================="
        require_env FFTW_LDFLAGS
        require_env FFTW_LIBS
        require_env FFTW_CFLAGS
        require_env ELPAROOT
        require_env ELPA_LDFLAGS
        require_env ELPA_LIBS
        require_env ELPA_CFLAGS
        require_env GSL_LDFLAGS
        require_env GSL_CFLAGS
        require_env GSL_LIBS
        require_env MATH_LIBS
        require_env MPI_LDFLAGS
        require_env MPI_LIBS
        require_env SCALAPACKROOT
        require_env SCALAPACK_LDFLAGS
        require_env SCALAPACK_CFLAGS
        require_env SCALAPACK_LIBS
        require_env LIBXC_LIBS
        require_env LIBXC_CFLAGS
        require_env LIBXC_LDFLAGS
        require_env SPGLIB_LIBS
        require_env SPGLIB_CFLAGS
        require_env SPGLIB_LDFLAGS
        require_env HDF5_LIBS
        require_env HDF5_CFLAGS
        require_env HDF5_LDFLAGS
        require_env LIBVDWXC_CFLAGS
        require_env LIBVDWXC_LIBS
        require_env LIBVDWXC_LDFLAGS
        require_env GSL_INCLUDE_DIR
        require_env GSL_LIBRARY
        require_env GSL_CBLAS_LIBRARY

        ARCH=`arch`
        SIRIUS_OPT="-O3 -DNDEBUG -mtune=native -ftree-loop-vectorize ${MATH_CFLAGS}"
        if [ "$ARCH" = "ppc64le" ]; then
            SIRIUS_OPT="-O3 -DNDEBUG -mcpu=power8 -mtune=power8 -funroll-loops -ftree-vectorize  -mvsx  -maltivec  -mpopcntd  -mveclibabi=mass -fvect-cost-model -fpeel-loops -mcmodel=medium ${MATH_CFLAGS}"
            SIRIUS_DBG="-O2 -g -mcpu=power8 -mtune=power8 -funroll-loops -ftree-vectorize  -mvsx  -maltivec  -mpopcntd  -mveclibabi=mass -fvect-cost-model -fpeel-loops -mcmodel=medium ${MATH_CFLAGS}"
        fi

        if [ "$ARCH" = "x86_64" ]; then
            SIRIUS_OPT="-O3 -DNDEBUG -mtune=native -ftree-loop-vectorize ${MATH_CFLAGS}"
            SIRIUS_DBG="-O2 -g -mtune=native -ftree-loop-vectorize ${MATH_CFLAGS}"
        fi

        pkg_install_dir="${INSTALLDIR}/sirius-${sirius_ver}"
        install_lock_file="${pkg_install_dir}/install_successful"
        if verify_checksums "${install_lock_file}" ; then
            echo "sirius_dist-${sirius_ver} is already installed, skipping it."
        else
            if [ -f SIRIUS-${sirius_ver}.tar.gz ] ; then
                echo "sirius_${sirius_ver}.tar.gz is found"
            else

                download_pkg ${DOWNLOADER_FLAGS} \
                             https://www.cp2k.org/static/downloads/SIRIUS-${sirius_ver}.tar.gz
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d sirius-${sirius_ver} ] && rm -rf sirius-${sirius_ver}
            tar -xzf SIRIUS-${sirius_ver}.tar.gz
            cd SIRIUS-${sirius_ver}
            cd src && sed -i -e "s/desciption/description/g" options.json && cd ..
            rm -Rf build
            mkdir build
            cd build
            COMPILATION_OPTIONS="-DHDF5_DIR=${HDF5_DIR}"
            if [ -n "$ELPA_LIBS" ] ; then
                if [ -s "$ELPAROOT" ] ; then
                    export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$ELPAROOT/lib/pkgconfig:$ELPAROOT/lib64/pkgconfig
                fi
                COMPILATION_OPTIONS="-DUSE_ELPA=ON -DELPA_INCLUDE_DIR=${ELPAROOT}/include/elpa-${ELPAVERSION}/elpa $COMPILATION_OPTIONS"
            fi

            if [ -n "$SCALAPACK_LIBS" ] ; then
                export SCALAPACK_LIB="$SCALAPACK_LIBS"
                if [ -s "$SCALAPACKROOT" ] ; then
                    COMPILATION_OPTIONS="-DUSE_SCALAPACK=ON -DSCALAPACK_INCLUDE_DIR=${SCALAPACKROOT}/include ${COMPILATION_OPTIONS}"
                else
                    COMPILATION_OPTIONS="-DUSE_SCALAPACK=ON ${COMPILATION_OPTIONS}"
                fi
            fi
            if [ -n "$HDF5_LIBS" ] ; then
                COMPILATION_OPTIONS="-DUSE_HDF5=ON $COMPILATION_OPTIONS"
            fi
            if [ -n "$LIBVDWXC_LIBS" ] ; then
                COMPILATION_OPTIONS="-DUSE_VDWXC=ON $COMPILATION_OPTIONS"
            fi
            if [ -n "$MKL_LIBS" ] ; then
                COMPILATION_OPTIONS="-DUSE_MKL=ON -DUSE_SCALAPACK=ON $COMPILATION_OPTIONS"
            fi

            echo "cmake -DCMAKE_INSTALL_PREFIX=${pkg_install_dir} \
                  -DCMAKE_CXXFLAGS_RELEASE="${SIRIUS_OPT}" \
                  -DCMAKE_CXX_FLAGS_RELWITHDEBINFO="${SIRIUS_DBG}" \
                  -DCMAKE_CXX_COMPILER=mpic++ \
                  -DCMAKE_C_COMPILER=mpicc ${COMPILATION_OPTIONS} .."

             cmake -DCMAKE_INSTALL_PREFIX=${pkg_install_dir} \
                  -DCMAKE_CXXFLAGS_RELEASE="${SIRIUS_OPT}" \
                  -DCMAKE_CXX_FLAGS_RELWITHDEBINFO="${SIRIUS_DBG}" \
                  -DCMAKE_CXX_COMPILER=mpic++ \
                  -DCMAKE_C_COMPILER=mpicc \
                  ${COMPILATION_OPTIONS} .. > compile.log 2>&1
             make -C src >> compile.log 2>&1

            install -d ${pkg_install_dir}/include >> install.log 2>&1
            install -d ${pkg_install_dir}/lib >> install.log 2>&1
            cp -R ../src/* ${pkg_install_dir}/include >> install.log 2>&1
            install -m 644 src/*.a ${pkg_install_dir}/lib >> install.log 2>&1
            install -m 644 src/mod_files/*.mod ${pkg_install_dir}/include >> install.log 2>&1
            cd ..

            # now do we have cuda as well

            if [ "$ENABLE_CUDA" = "__TRUE__" ] ; then
                [ -d build-cuda ] && rm -rf "build-cuda"
                mkdir build-cuda
                cd build-cuda
                cmake -DCMAKE_INSTALL_PREFIX=${pkg_install_dir} \
                      -DCMAKE_CXXFLAGS_RELEASE="${SIRIUS_OPT}" \
                      -DCMAKE_CXX_FLAGS_RELWITHDEBINFO="${SIRIUS_DBG}" \
                      -DUSE_CUDA=ON \
                      -DGPU_MODEL=P100 \
                      -DCMAKE_CXX_COMPILER=mpic++ \
                      -DCMAKE_C_COMPILER=mpicc ${COMPILATION_OPTIONS} .. >> compile.log 2>&1
                make -C src >> compile.log 2>&1
                install -m 644 src/*.a ${pkg_install_dir}/lib/cuda >> install.log 2>&1
                install -m 644 src/mod_files/*.mod ${pkg_install_dir}/include/cuda >> install.log 2>&1
                SIRIUS_CUDA_LDFLAGS="-L'${pkg_install_dir}/lib/cuda' -Wl,-rpath='${pkg_install_dir}/lib/cuda'"
                cd ..
            fi
            SIRIUS_CFLAGS="-I'${pkg_install_dir}/include/cuda'"
            SIRIUS_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
            write_checksums "${install_lock_file}" "${SCRIPT_DIR}/$(basename ${SCRIPT_NAME})"
        fi
        ;;
    __SYSTEM__)
        require_env FFTW_LDFLAGS
        require_env FFTW_LIBS
        require_env FFTW_CFLAGS
        require_env ELPAROOT
        require_env ELPA_LDFLAGS
        require_env ELPA_LIBS
        require_env ELPA_CFLAGS
        require_env GSL_LDFLAGS
        require_env GSL_CFLAGS
        require_env GSL_LIBS
        require_env MATH_LIBS
        require_env MPI_LDFLAGS
        require_env MPI_LIBS
        require_env SCALAPACKROOT
        require_env SCALAPACK_LDFLAGS
        require_env SCALAPACK_CFLAGS
        require_env SCALAPACK_LIBS
        require_env LIBXC_LIBS
        require_env LIBXC_CFLAGS
        require_env LIBXC_LDFLAGS
        require_env SPGLIB_LIBS
        require_env SPGLIB_CFLAGS
        require_env SPGLIB_LDFLAGS
        require_env HDF5_LIBS
        require_env HDF5_CFLAGS
        require_env HDF5_LDFLAGS
        require_env LIBVDWXC_CFLAGS
        require_env LIBVDWXC_LDFLAGS
        require_env LIBVDWXC_LIBS

        check_lib -lsirius_f "sirius_f"
        add_include_from_paths SIRIUS_CFLAGS "sirius*" $INCLUDE_PATHS
        add_lib_from_paths SIRIUS_LDFLAGS "libsirius_f.*" $LIB_PATHS
        ;;
    *)
        echo "==================== Linking SIRIUS_Dist to user paths ===================="
        pkg_install_dir="$with_sirius"
        check_dir "${pkg_install_dir}/lib"
        check_dir "${pkg_install_dir}/lib64"
        check_dir "${pkg_install_dir}/include"
        ;;
esac
if [ "$with_sirius" != "__DONTUSE__" ] ; then
    SIRIUS_LIBS="-lsirius_f"
    SIRIUS_CUDA_LDFLAGS="-L'${pkg_install_dir}/lib/cuda' -Wl,-rpath='${pkg_install_dir}/lib/cuda'"
    SIRIUS_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
    SIRIUS_CFLAGS="-I'${pkg_install_dir}/include'"
    if [ "$with_sirius" != "__SYSTEM__" ] ; then
        cat <<EOF > "${BUILDDIR}/setup_sirius"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib/cuda"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib/cuda"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib/cuda"
prepend_path CPATH "$pkg_install_dir/include"
EOF
        cat "${BUILDDIR}/setup_sirius" >> $SETUPFILE
    fi
    cat <<EOF >> "${BUILDDIR}/setup_sirius"
export SIRIUS_CFLAGS="IF_CUDA(-I${pkg_install_dir}/include/cuda|-I${pkg_install_dir}/include)"
export SIRIUS_FFLAGS="IF_CUDA(-I${pkg_install_dir}/include/cuda|-I${pkg_install_dir}/include)"
export SIRIUS_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
export SIRIUS_CUDA_LDFLAGS="-L'${pkg_install_dir}/lib/cuda' -Wl,-rpath='${pkg_install_dir}/lib/cuda'"
export SIRIUS_LIBS="${SIRIUS_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(IF_OMP("-D__SIRIUS"|)|)"
export CP_CFLAGS="\${CP_CFLAGS} IF_MPI(IF_OMP("\${SIRIUS_CFLAGS}"|)|)"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_MPI(IF_OMP(IF_CUDA("\${SIRIUS_CUDA_LDFLAGS}"|"\${SIRIUS_LDFLAGS}")|)|)"
export CP_LIBS="IF_MPI(IF_OMP("\${SIRIUS_LIBS}"|)|) \${CP_LIBS}"
EOF
fi
cd "${ROOTDIR}"
source "${SCRIPT_DIR}"/step_end.sh