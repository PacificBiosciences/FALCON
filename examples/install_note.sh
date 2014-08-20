# This is the script that will build everything needed to generate an assembly 
# on top of the StarCluster Ubuntu AMI 

HBAR_ROOT=$HOME #if you don't want to install this in your home directory, replace $HOME with a proper directory
mkdir -p $HBAR_ROOT/HBAR_ENV
export HBAR_HOME=$HBAR_ROOT/HBAR_ENV/
sudo apt-get install python-virtualenv
virtualenv -p /usr/bin/python2.7 $HBAR_HOME
cd $HBAR_HOME
. bin/activate
pip install numpy==1.6.2
sudo apt-get install python-dev
pip install numpy==1.6.2
wget http://www.hdfgroup.org/ftp/HDF5/prev-releases/hdf5-1.8.9/src/hdf5-1.8.9.tar.gz
tar zxvf hdf5-1.8.9.tar.gz
cd hdf5-1.8.9
./configure --prefix=$HBAR_HOME --enable-cxx
make install
cd ..
wget http://h5py.googlecode.com/files/h5py-2.0.1.tar.gz
tar zxvf h5py-2.0.1.tar.gz
cd h5py-2.0.1
python setup.py build --hdf5=$HBAR_HOME
python setup.py install
cd ..
pip install git+https://github.com/PacificBiosciences/pbcore.git#pbcore
sudo apt-get install git
pip install git+https://github.com/PacificBiosciences/pbcore.git#pbcore
pip install git+https://github.com/PacificBiosciences/pbdagcon.git#pbdagcon
pip install git+https://github.com/PacificBiosciences/pbh5tools.git#pbh5tools
pip install git+https://github.com/cschin/pypeFLOW.git#pypeflow
pip install rdflib==3.4.0
pip install git+https://github.com/PacificBiosciences/HBAR-DTK.git#hbar-dtk
pip install git+https://github.com/PacificBiosciences/FALCON.git#falcon

git clone https://github.com/PacificBiosciences/blasr.git
cd blasr
export HDF5INCLUDEDIR=/home/HBAR_ENV/include/
export HDF5LIBDIR=/home/HBAR_ENV/lib/
make
cp alignment/bin/blasr ../bin/
cp alignment/bin/sawriter ../bin/
cp pbihdfutils/bin/samFilter  ../bin
cp pbihdfutils/bin/samtoh5  ../bin
cd ..


wget http://downloads.sourceforge.net/project/boost/boost/1.47.0/boost_1_47_0.tar.gz
tar zxvf boost_1_47_0.tar.gz
cd boost_1_47_0/
bash bootstrap.sh
./b2 install -j 24 --prefix=$HBAR_ROOT/HBAR_ENV/boost
cd ..

sudo apt-get install libpcre3 libpcre3-dev
wget http://downloads.sourceforge.net/project/swig/swig/swig-2.0.11/swig-2.0.11.tar.gz
tar zxvf swig-2.0.11.tar.gz
cd swig-2.0.11
./configure --prefix=$HBAR_ROOT/HBAR_ENV
make
make install
cd ..

git clone https://github.com/PacificBiosciences/ConsensusCore.git
cd ConsensusCore/
python setup.py install --swig=$HBAR_ROOT/HBAR_ENV/bin/swig --boost=$HBAR_ROOT/HBAR_ENV/boost/include/
cd ..

pip install git+https://github.com/PacificBiosciences/GenomicConsensus.git#GenomicConsensus
pip install git+https://github.com/PacificBiosciences/pbalign#pbaligno

wget http://downloads.sourceforge.net/project/mummer/mummer/3.23/MUMmer3.23.tar.gz
tar zxvf MUMmer3.23.tar.gz
cd MUMmer3.23/
make install
cd ..
export PATH=$PATH:/home/HBAR_ENV/MUMmer3.23


wget http://downloads.sourceforge.net/project/samtools/samtools/0.1.19/samtools-0.1.19.tar.bz2
tar jxvf samtools-0.1.19.tar.bz2
cd samtools-0.1.19
make
cp samtools ../bin
cd ..
