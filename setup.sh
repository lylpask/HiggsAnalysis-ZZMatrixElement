#!/bin/sh


#echo "$@"
pushd ${CMSSW_BASE}/src/ZZMatrixElement/MELA/
source setup.sh "$@"
popd
scramv1 b "$@"
