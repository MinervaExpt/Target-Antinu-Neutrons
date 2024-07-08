#!/bin/bash
if [ -n "$ROOTSYS" ]; then
  echo "ROOT is already set up.  You must set up the MINERvA 101 tutorial from scratch.  Check your .bash_profile for \"root\""
fi

source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
spack load root@6.28.12
spack load cmake
spack load gcc
#spack load fife-utils Not sure if needed

export LD_LIBRARY_PATH=${ROOTSYS}/lib/root:${LD_LIBRARY_PATH}
