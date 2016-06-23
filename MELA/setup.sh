#!/bin/sh

#echo "$@"
if [[ "$1" == *"clean"* ]];then
	scramv1 b "$@"
	pushd ${CMSSW_BASE}/src/ZZMatrixElement/MELA/fortran/
	make clean
	rm -f ../data/${SCRAM_ARCH}/libjhugenmela.so
	popd
else
	pushd ${CMSSW_BASE}/src/ZZMatrixElement/MELA/fortran/
	make all
	mv libjhugenmela.so ../data/${SCRAM_ARCH}/
	popd
	scramv1 b "$@"
fi
