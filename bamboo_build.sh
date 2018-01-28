#!/bin/bash

type module >& /dev/null || . /mnt/software/Modules/current/init/bash
module purge
module load git
module load gcc
module load ccache
module load python/2
#module load make

set -vx
git --version
which gcc
which g++
gcc --version
which python

mkdir -p LOCAL
export PYTHONUSERBASE=$(pwd)/LOCAL
export PATH=${PYTHONUSERBASE}/bin:${PATH}

# We need latest local pypeFLOW, not from PyPI.
pushd ../pypeFLOW
pop install --user --edit .
popd
# Back to FALCON.

make install-edit
# Note: no --edit because we might be building artifacts.
# ... Scratch that. We have trouble getting coverage for
#  source=falcon_kit
# but maybe it will work with a --edit install.

export MY_TEST_FLAGS="-v -s --durations=0 --cov=falcon_kit --cov-report=term-missing --cov-report=xml:coverage.xml --cov-branch"
make test
sed -i -e 's@filename="@filename="./falcon_kit/@g' coverage.xml

make pylint
