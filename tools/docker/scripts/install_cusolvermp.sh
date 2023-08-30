#!/bin/bash -e

# author: Ole Schuett

if (($# != 1)); then
  echo "Usage: install_cusolvermp.sh <GPU_VERSION>"
  exit 1
fi

GPU_VERSION=$1

case ${GPU_VERSION} in
  P100)
    NVCC_GENCODE="--with-nvcc-gencode='-arch sm_60'"
    ;;
  V100)
    NVCC_GENCODE="--with-nvcc-gencode='-arch sm_70'"
    ;;
  A100)
    NVCC_GENCODE="--with-nvcc-gencode='-arch sm_80'"
    ;;
  *)
    NVCC_GENCODE="" # Build CUDA kernels for all GPU versions, which takes a long time.
    ;;
esac

# Install UCX
echo -e "Installing UCX..."
wget -q https://github.com/openucx/ucx/releases/download/v1.14.1/ucx-1.14.1-ubuntu22.04-mofed5-cuda11.tar.bz2
echo "d5f7f6e46a1495edf232cee63ba71fa3f60bfb31f8cd69ae7b62e84f1ea01649  ucx-1.14.1-ubuntu22.04-mofed5-cuda11.tar.bz2" | sha256sum --check
tar -xjf ucx-1.14.1-ubuntu22.04-mofed5-cuda11.tar.bz2
dpkg -i ucx-1.14.1.deb ucx-cuda-1.14.1.deb
rm -rf ucx-*

# Install UCC
echo -e "Installing UCC..."
wget -q https://github.com/openucx/ucc/archive/refs/tags/v1.2.0.tar.gz -O ucc-1.2.0.tar.gz
echo "c1552797600835c0cf401b82dc89c4d27d5717f4fb805d41daca8e19f65e509d  ucc-1.2.0.tar.gz" | sha256sum --check
tar -xzf ucc-1.2.0.tar.gz
cd ucc-1.2.0
./autogen.sh &> autogen.log
./configure --prefix=/usr --with-cuda="$CUDA_PATH" --with-ucx=/usr "${NVCC_GENCODE}" &> config.log
make -j20 &> build.log
make install &> install.log
cd ..
rm -rf ucc-*

# Install cuSolverMp
echo -e "Installing cuSolverMp..."
wget -q https://developer.download.nvidia.com/assets/gameworks/downloads/regular/cuSolverMp/cusolvermp-linux-x86_64-v0.4.1.0-cuda11.8.tar.gz
echo "69f374d36578e95445a877eb928a1a98942b4323680aa343a9cd50a2611c5aec  cusolvermp-linux-x86_64-v0.4.1.0-cuda11.8.tar.gz" | sha256sum --check
tar -xzf cusolvermp-linux-x86_64-v0.4.1.0-cuda11.8.tar.gz
cd cusolvermp-linux-x86_64-v0.4.1.0-cuda11.8
cp -r include/* /usr/include/
cp -r lib/* /usr/lib/
#TODO Find a way to move ucc.conf from /share to /usr/share.
mkdir /share
cp share/ucc.conf /share/
cd ..
rm -r cusolvermp-*

# EOF
