#!/bin/sh

export PATH=/mnt/software/a/anaconda2/4.2.0/bin:$PATH
which python

mkdir -p LOCAL
export PYTHONUSERBASE=$(pwd)/LOCAL

pip install --user nose coverage
nosetests -v --with-xunit --xunit-file=nose.basic.xml --with-coverage --cover-xml --cover-xml-file=coverage.xml test/
sed -i -e 's@filename="@filename="./@g' coverage.xml
nosetests -v --with-xunit --xunit-file=nose.doctest.xml --with-doctest falcon_kit/functional.py
