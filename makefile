# Feel free to override this.
ifndef PYTHONUSERBASE
  PYTHONUSERBASE:=LOCAL
  PATH:=${PYTHONUSERBASE}/bin:${PATH}
  export PYTHONUSERBASE
  export PATH
endif
MY_TEST_FLAGS?=-v -s

install-edit:
	pip -v install --user --edit .
install: wheel
	pip -v install --user --use-wheel --find-links=dist/ .
test:
	python -c 'import falcon_kit; print falcon_kit.falcon'
	#pip install --user pytest coverage
	coverage run --source falcon_kit -m py.test ${MY_TEST_FLAGS} --junit-xml=test.basic.xml test/
	coverage run --source falcon_kit -m py.test ${MY_TEST_FLAGS} --junit-xml=test.doctest.xml --doctest-modules falcon_kit/functional.py
	cp -f test.basic.xml nose.basic.xml
	cp -f test.doctest.xml nose.doctest.xml
	coverage xml -o coverage.xml
	sed -i -e 's@filename="@filename="./@g' coverage.xml

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
