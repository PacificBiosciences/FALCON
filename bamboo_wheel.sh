#!/bin/bash
type module >& /dev/null || . /mnt/software/Modules/current/init/bash
set -vex
ls -larth ..
ls -larth
pwd

module purge
module load gcc
module load ccache

module load python/2-UCS2
make wheel WHEELHOUSE=artifacts/gcc-6.4.0/wheelhouse

# For now, we have only "any" wheels, so we do not need to build again.

module unload python

module load python/2-UCS4
make wheel WHEELHOUSE=artifacts/gcc-6.4.0/wheelhouse
