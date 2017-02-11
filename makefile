# Feel free to override this.
ifndef PYTHONUSERBASE
  PYTHONUSERBASE:=LOCAL
  PATH:=${PYTHONUSERBASE}/bin:${PATH}
  export PYTHONUSERBASE
  export PATH
endif

install-edit:
	pip -v install --user --edit .
install: wheel
	pip -v install --user --use-wheel --find-links=dist/ .
test:
	python -c 'import falcon_kit; print falcon_kit.falcon'
	pip install --user nose
	nosetests -v --with-xunit --xunit-file=nose.basic.xml test/
	nosetests -v --with-xunit --xunit-file=nose.doctest.xml --with-doctest falcon_kit/functional.py
# We cannot run doctests on *all* modules because some include dependencies.
# Just pypeFLOW for now, but I would rather not test dependencies anyway.

wheel:
	pip install --upgrade --user pip
	python setup.py bdist_wheel
# Look for dist/*.whl

tar:
	rm -f FALCON.tar.gz
	tar cvzf FALCON.tar.gz -C ${PYTHONUSERBASE} .
# Much smaller than the wheel, and includes all necessary dependencies,
# but also includes anything already in the user-site.


.PHONY: install test install-no-edit wheel
