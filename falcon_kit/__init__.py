from __future__ import absolute_import

from .falcon_kit import *

try:
    import sys, pkg_resources
    sys.stderr.write('{}\n'.format(pkg_resources.get_distribution('falcon-kit')))
except Exception:
    pass
