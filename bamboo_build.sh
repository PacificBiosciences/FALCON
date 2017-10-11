#!/bin/bash

type module >& /dev/null || . /mnt/software/Modules/current/init/bash
module unload git gcc ccache
module load git/2.8.3
module load gcc/6.4.0
module load ccache/3.2.3
#module load make

set -vx
git --version
which gcc
which g++
gcc --version
# We cannot use /bin/python without /bin/gcc.
export PATH=/mnt/software/a/anaconda2/4.2.0/bin:$PATH
which python

mkdir -p LOCAL
export PYTHONUSERBASE=$(pwd)/LOCAL

# We have a problem with pylint: https://github.com/PyCQA/pylint/issues/1296
pip install --user --upgrade pylint

make install-edit
# Note: no --edit because we might be building artifacts.
# ... Scratch that. We have trouble getting coverage for
#  source=falcon_kit
# but maybe it will work with a --edit install.

make pylint
