#!/bin/sh
set -vex
export PATH=/mnt/software/a/anaconda2/4.2.0/bin:$PATH
which python

mkdir -p LOCAL
export PYTHONUSERBASE=$(pwd)/LOCAL
export PATH=$PYTHONUSERBASE/bin:$PATH

pip install --user pytest coverage
#make test
make coverage-install
make coverage
