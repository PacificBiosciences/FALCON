#!/bin/sh
# -e: fail on error
# -v: show commands
# -x: show expanded commands
set -vex

#env | sort
mkdir -p LOCAL
export PYTHONUSERBASE=$(pwd)/LOCAL
export PATH=$PYTHONUSERBASE/bin:$PATH

git clone git://github.com/PacificBiosciences/pypeFLOW.git
cd pypeFLOW
git checkout origin/develop
pip install --user --edit .
cd ..
make install
make test
