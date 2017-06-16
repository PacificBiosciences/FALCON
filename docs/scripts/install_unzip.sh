#!/bin/bash
# This is a simple script to install FALCON-integrate + FALCON_unzip.
# Running FALCON_unzip requires addiitonal blasr/samtools/GenomicConsensus binaries
# contact: gconcepcion@pacificbiosciences.com

###Install script dependencies
##This script should work on both Ubuntu/CentOS as long as the following dependencies are installed and 
##available in your $PATH

#source /mnt/software/Modules/current/init/bash
#module load python/2.7.9
#module load virtualenv/13.0.1
#module load git


###Variables
ROOT=$1
if [ -z $ROOT ];then
    ROOT=$(pwd)
fi

if [ ! -z $ROOT ]; then
    ROOT=$(readlink -f $ROOT)
fi

SRC_ROOT=${ROOT}/src
VENV_BASE=${ROOT}/fc_env
SLUG=$(date +%y%m%d)
VENV=${VENV_BASE}_${SLUG}
REPO_ROOT="http://github.com/PacificBiosciences"
UNZIP="FALCON_unzip"
FALCON="FALCON-integrate"
UNZIP_REPO=${REPO_ROOT}/${UNZIP}.git
FALCON_REPO=${REPO_ROOT}/${FALCON}.git
UNZIP_PATH=${SRC_ROOT}/${UNZIP}
FALCON_PATH=${SRC_ROOT}/${FALCON}

###Cleanup
if [ ! -d "${SRC_ROOT}" ]; then
    echo "Creating src root"
    mkdir ${SRC_ROOT}
fi

if [ -d "${VENV}" ]; then
    echo "Removing previous build from today"
    rm -rf $VENV
fi

if [ -d "${UNZIP_PATH}" ]; then
    echo "Removing old UNZIP repo"
    rm -rf ${UNZIP_PATH}
fi

if [ -d "${FALCON_PATH}" ]; then
    echo "Removing old FALCON-integrate repo"
    rm -rf ${FALCON_PATH}
fi

if [ -L "${VENV_BASE}" ]; then
    unlink ${VENV_BASE}
fi


cd $SRC_ROOT
git clone ${FALCON_REPO}
git clone ${UNZIP_REPO}


###install falcon
virtualenv --no-site-packages ${VENV}
source ${VENV}/bin/activate
cd ${FALCON_PATH}
git checkout master
git submodule update --init
sed -i "s|^FALCON_WORKSPACE.*|FALCON_WORKSPACE=${FALCON_PATH}|g" default-env.sh
sed -i "s|^PYTHONUSERBASE.*|PYTHONUSERBASE=${VENV}|g" default-env.sh
make init
source env.sh

make config-standard
make -j all
make install
#make test


###install unzip
cd ${UNZIP_PATH}
python setup.py install
ln -s ${VENV} ${VENV_BASE}


###Test falcon_unzip pipeline
## requires blasr, samtools and GenomicConsensus binaries in your path
#cd ${FALCON_PATH}/FALCON-examples
#../git-sym/git-sym update run/greg200k-sv2
#cd run/greg200k-sv2
#fc_run fc_run.cfg
#fc_unzip.py fc_unzip.cfg
