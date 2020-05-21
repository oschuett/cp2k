#!/bin/bash -e

# author: Ole Schuett

# Install Ubuntu packages required for the toolchain.

echo "Installing Ubuntu packages..."

export DEBIAN_FRONTEND=noninteractive
export DEBCONF_NONINTERACTIVE_SEEN=true

apt-get update -qq

apt-get install -qq --no-install-recommends \
    autoconf                                \
    autogen                                 \
    automake                                \
    autotools-dev                           \
    ca-certificates                         \
    g++                                     \
    git                                     \
    less                                    \
    libtool                                 \
    make                                    \
    nano                                    \
    pkg-config                              \
    python3-numpy                           \
    python3                                 \
    unzip                                   \
    wget                                    \
    xxd                                     \
    zlib1g-dev

rm -rf /var/lib/apt/lists/*

#EOF
