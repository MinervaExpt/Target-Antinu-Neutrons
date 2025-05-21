#!/bin/bash
if [ -n "$ROOTSYS" ]; then
  echo "ROOT is already set up.  You must set up the MINERvA 101 tutorial from scratch.  Check your .bash_profile for \"root\""
fi

source /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh
spack load root@6.28.12 arch=linux-almalinux9-x86_64_v3
spack load cmake arch=linux-almalinux9-x86_64_v3
spack load gcc
#spack load fife-utils Not sure if needed

export LD_LIBRARY_PATH=${ROOTSYS}/lib/root:${LD_LIBRARY_PATH}
