virtualenv --no-site-packages  --always-copy   $PWD/fc_env
. $PWD/fc_env/bin/activate

cd FALCON
git submodule init
git submodule update

cd pypeFLOW
python setup.py install
cd ..

python setup.py install

cd DAZZ_DB/
make
cp DBrm DBshow DBsplit DBstats fasta2DB ../../fc_env/bin/
cd ..

cd DALIGNER
make
cp daligner daligner_p DB2Falcon HPCdaligner LA4Falcon LAmerge LAsort  ../../fc_env/bin
cd ../..
