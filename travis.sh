#!/bin/sh
# -e: fail on error
# -v: show commands
# -x: show expanded commands
set -vex

#env | sort
mkdir -p LOCAL
export PYTHONUSERBASE=$(pwd)/LOCAL
export PATH=$PYTHONUSERBASE/bin:$PATH

. ./install.sh
. ./test.sh
