#!/bin/bash

if [ -z "${BASH_SOURCE[0]}" ]; then
  export DRICH_F4A=$(dirname $(realpath $0))
else
  export DRICH_F4A=$(dirname $(realpath ${BASH_SOURCE[0]}))
fi
echo "DRICH_F4A=$DRICH_F4A"
export DRICH_HOME=$DRICH_F4A/..

dirInstall=$DRICH_F4A/install
if [ -d "$dirInstall" ]; then
  source /cvmfs/eic.opensciencegrid.org/ecce/gcc-8.3/opt/fun4all/core/bin/setup_local.sh $dirInstall
fi
