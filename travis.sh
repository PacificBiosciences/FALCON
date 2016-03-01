#!/bin/sh
# -e: fail on error
# -v: show commands
# -x: show expanded commands
set -vex

#env | sort
mkdir -p fc-env
rm -f fc-env/bin/python
virtualenv -p python2.7 fc-env || ../virtualenv/virtualenv.py fc-env
. fc-env/bin/activate
python setup.py -v install
python -c 'import falcon_kit; print falcon_kit.falcon'

# When doctests are passing, add this:
pip install nose
nosetests -v test/
nosetests -v --with-doctest falcon_kit/functional.py
# We cannot run that on *all* modules because some include dependencies.
# Just pypeFLOW for now, but I would rather not test dependencies.
