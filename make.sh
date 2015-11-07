#!/bin/bash
source ${RESCOL_DIR}/common/variable.sh
build_dir=${RESCOL_DIR}/build/${PETSC_ARCH}

if [ ! -e $build_dir ]; then
  mkdir -p $build_dir
fi
cd $build_dir
pwd
make --file=${RESCOL_DIR}/makefile $@
