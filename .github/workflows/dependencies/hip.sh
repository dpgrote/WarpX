#!/usr/bin/env bash
#
# Copyright 2020 The WarpX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

# Ref.: https://rocmdocs.amd.com/en/latest/Installation_Guide/Installation-Guide.html#ubuntu
wget -q -O - http://repo.radeon.com/rocm/rocm.gpg.key \
  | sudo apt-key add -
echo 'deb [arch=amd64] http://repo.radeon.com/rocm/apt/debian/ xenial main' \
  | sudo tee /etc/apt/sources.list.d/rocm.list

echo 'export PATH=$PATH:/opt/rocm/bin:/opt/rocm/profiler/bin:/opt/rocm/opencl/bin' \
  | sudo tee -a /etc/profile.d/rocm.sh
# we should not need to export HIP_PATH=/opt/rocm/hip with those installs

sudo apt-get update

# Ref.: https://rocmdocs.amd.com/en/latest/Installation_Guide/Installation-Guide.html#installing-development-packages-for-cross-compilation
# meta-package: rocm-dkms
# OpenCL: rocm-opencl
# other: rocm-dev rocm-utils
sudo apt-get install -y --no-install-recommends \
    build-essential \
    gfortran        \
    libnuma-dev     \
    libopenmpi-dev  \
    openmpi-bin     \
    rocm-dev        \
    rocfft          \
    rocrand

# activate
#
source /etc/profile.d/rocm.sh
hipcc --version

# cmake-easyinstall
#
sudo curl -L -o /usr/local/bin/cmake-easyinstall https://git.io/JvLxY
sudo chmod a+x /usr/local/bin/cmake-easyinstall
export CEI_SUDO="sudo"
