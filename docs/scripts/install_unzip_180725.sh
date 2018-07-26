#!/bin/bash
# This is a simple script to install FALCON-integrate + FALCON_unzip.
# contact: gconcepcion@pacificbiosciences.com

###Install script dependencies
##This script should work on both Ubuntu/CentOS as long as the following dependencies are installed and
##available in your $PATH

source /mnt/software/Modules/current/init/bash

#module load python/2.7.14-UCS4
#module load virtualenv/13.0.1
###These three are necessary for Unzip
#module load samtools
#module load minimap2
#module load gcc/6.4.0
#module load git/2.8.3

###Variables
INPUT=$1
if [ -z ${INPUT} ]; then
    ROOT=$(pwd)
else
    ROOT=$(readlink -f ${INPUT})
fi

SRC=${ROOT}/src
VENV_BASE=${ROOT}/fc_env
SLUG=$(date +%y%m%d)
VENV=${VENV_BASE}_${SLUG}
FALCON_PATH=${SRC}/FALCON-integrate
FALCON_REPO="https://github.com/PacificBiosciences/FALCON-integrate.git"
DL_CLOUD="https://downloads.pacbcloud.com/public/falcon"

python --version
###Test whether you need ucs4 or ucs2 tarball with this command:
UCS_VER=$(python2.7 -c 'import sysconfig,pprint; pprint.pprint(sysconfig.get_config_vars()["Py_UNICODE_SIZE"])')
UNZIP_TARBALL="falcon-2018.07.25-07.01-py2.7-ucs${UCS_VER}-beta.tar.gz"
###Cleanup
if [ ! -d ${ROOT} ]; then
    mkdir ${ROOT}
fi

if [ -d ${SRC} ]; then
    rm -rf ${SRC}
    mkdir ${SRC}
else
    mkdir ${SRC}
fi 

if [ -d "${VENV}" ]; then
    echo "Removing previous build from today"
    rm -rf $VENV
fi

if [ -d "${FALCON_PATH}" ]; then
    echo "Removing old FALCON-integrate repo"
    rm -rf ${FALCON_PATH}
fi

if [ -L "${VENV_BASE}" ]; then
    unlink ${VENV_BASE}
fi

###Make python virtualenv
virtualenv --python=$(which python) --no-site-packages ${VENV}
echo "export LD_LIBRARY_PATH=${VENV}/lib:\${LD_LIBRARY_PATH}" >>${VENV}/bin/activate
source ${VENV}/bin/activate

###Install FALCON and FALCON_unzip

cd ${SRC}
curl -v -k -O ${DL_CLOUD}/${UNZIP_TARBALL}
tar zxvf ${UNZIP_TARBALL} -C ${VENV}

#ln -s ${VENV} ${VENV_BASE}

###Install MUMmer 3.23 - for FALCON_unzip

MUMMER_323='https://downloads.sourceforge.net/project/mummer/mummer/3.23/MUMmer3.23.tar.gz'
cd ${SRC}
wget ${MUMMER_323} -P ${SRC} --no-check-certificate
tar zxvf MUMmer3.23.tar.gz
cd MUMmer3.23
make install
find . -maxdepth 1 -perm -111 -type f -exec cp {} ${VENV}/bin \;

###Test falcon_unzip pipeline
cd ${SRC} && git clone ${FALCON_REPO}
cd ${FALCON_PATH} && git submodule update --init
cd ${FALCON_PATH}/FALCON-examples && ../git-sym/git-sym update run/greg200k-sv2
git checkout 7aa47ddd99a6c2c8a87e76bce83eaa320b3d9474
cd run/greg200k-sv2 && fc_run fc_run.cfg
fc_unzip.py fc_unzip.cfg

echo "FALCON & FALCON_unzip have been successfully installed into a virtualenv."
echo "To activate the FALCON_unzip environment:"
echo "$ source ${VENV}/bin/activate"
