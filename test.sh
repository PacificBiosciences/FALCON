pip install --user nose

nosetests -v test/

nosetests -v --with-doctest falcon_kit/functional.py
# We cannot run that on *all* modules because some include dependencies.
# Just pypeFLOW for now, but I would rather not test dependencies anyway.
