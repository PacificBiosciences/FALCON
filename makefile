# Feel free to override this.
ifndef PYTHONUSERBASE
  PYTHONUSERBASE:=LOCAL
  PATH:=${PYTHONUSERBASE}/bin:${PATH}
  export PYTHONUSERBASE
  export PATH
endif
export COVERAGE_PROCESS_START

MY_TEST_FLAGS?=-v -s

install-edit:
	pip -v install --user --edit .
install: wheel
	pip -v install --user --use-wheel --find-links=dist/ .
pylint:
	pylint --errors-only falcon_kit/
test:
	python -c 'import falcon_kit; print falcon_kit.falcon'
	#pip install --user pytest
	py.test ${MY_TEST_FLAGS} --junit-xml=test.basic.xml test/
	py.test ${MY_TEST_FLAGS} --junit-xml=test.doctest.xml --doctest-modules falcon_kit/functional.py
	cp -f test.basic.xml nose.basic.xml
	cp -f test.doctest.xml nose.doctest.xml
coverage:
	make coverage-clean
	#pip install --user coverage
	COVERAGE_PROCESS_START=${PWD}/mycoverage.cfg ${MAKE} coverage-actual
coverage-actual: test
	ls -larth
	coverage combine
	ls -larth
	coverage xml -o coverage.xml
	sed -i -e 's@filename="@filename="./@g' coverage.xml
	coverage report -m
coverage-clean:
	rm -f .coverage* coverage.xml
coverage-install:
	# This is needed only if you run from a different directory, since ./sitecustomize.py
	# would not be in 'sys.path'.
	# Assume PYTHONUSERBASE is set.
	mkdir -p ${PYTHONUSERBASE}/lib/python2.7/site-packages
	ln -f mysitecustomize.py ${PYTHONUSERBASE}/lib/python2.7/site-packages/sitecustomize.py
coverage-uninstall:
	rm -f ${PYTHONUSERBASE}/lib/python2.7/site-packages/sitecustomize.py*

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

clean: coverage-clean
	\rm -f *.xml


.PHONY: install test install-no-edit wheel coverage tar clean
