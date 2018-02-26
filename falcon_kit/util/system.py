from __future__ import absolute_import
from __future__ import unicode_literals

from future.utils import viewitems
from .io import system
import contextlib
import logging
import os
import pprint
import fnmatch

log = logging.getLogger(__name__)


def only_these_symlinks(dir2paths):
    """Create symlinks, and delete all other symlinks for each directory.
      dir2paths := {dir: [paths]}
    ('paths' is usually a list of 1.)
    Use relative symlink targets.
    Possibly, log everything we do, as one log statement to avoid user complaints.
    """
    log.info('Symlink .las files for further merging:\n{}'.format(
        pprint.pformat(dict(dir2paths))))
    for (d, paths) in viewitems(dir2paths):
        bases = [os.path.basename(path) for path in paths]
        base2rel = {os.path.basename(path): os.path.relpath(
            path, d) for path in paths}
        assert len(base2rel) == len(
            bases), 'Non-unique basename in {}'.format(repr(paths))
        for existing_base in os.listdir(d):
            existing_path = os.path.join(d, existing_base)
            if os.path.islink(existing_path):
                if existing_base in base2rel:
                    if os.readlink(existing_path) != base2rel[existing_base]:
                        # Wrong target (or non-relative) so remove it.
                        os.unlink(existing_path)
                    else:
                        del base2rel[existing_base]  # Just keep it.
                else:
                    os.unlink(existing_path)  # Old? Remove it for safety.
        for (base, rel) in viewitems(base2rel):
            path = os.path.join(d, base)
            os.symlink(rel, path)


def lfs_setstripe_maybe(path='.', stripe=12):
    path = os.path.abspath(path)
    rc = system('lfs setstripe -c {:d} {!s}'.format(stripe, path))
    if rc:
        log.info('Apparently {!r} is not lustre in filesystem.'.format(path))
    else:
        log.info('This lfs stripe ({}) should propagate to subdirs of {!r}.'.format(
            stripe, path))


def find_files(root_path, pattern):
    """
    Finds all files with filenames formatted as the
    given pattern, descending down from root_path.
    raise Exception if root_path is not a directory.
    """
    if not os.path.isdir(root_path):
        raise Exception('Not a directory: {!r}'.format(root_path))
    for root, dirs, files in os.walk(root_path):
        dirs.sort()
        for filename in sorted(fnmatch.filter(files, pattern)):
            yield os.path.join(root, filename)


@contextlib.contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    log.warning('CD: %r <- %r' % (newdir, prevdir))
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        log.warning('CD: %r -> %r' % (newdir, prevdir))
        os.chdir(prevdir)


def abs_fns(ifofns, idir=None):
    """Yield absolute filenames from a streamed file-of-filenames.
    """
    log.info('Absolutizing FOFN in dir={!r}'.format(os.path.abspath(idir)))
    for line in ifofns.read().split():
        ifn = line.strip()
        if not ifn:
            continue
        if not os.path.isabs(ifn):
            ifn = os.path.abspath(os.path.join(idir, ifn))
        yield ifn


def make_fofn_abs(i_fofn_fn, o_fofn_fn):
    """Copy i_fofn to o_fofn, but with relative filenames expanded for the dir of i_fofn.
    """
    assert os.path.abspath(o_fofn_fn) != os.path.abspath(
        i_fofn_fn), '{!r} != {!r}'.format(o_fofn_fn, i_fofn_fn)
    with open(i_fofn_fn) as ifs, open(o_fofn_fn, 'w') as ofs:
        for fn in abs_fns(ifs, os.path.dirname(os.path.realpath(i_fofn_fn))):
            ofs.write(fn + '\n')
    # return o_fofn_fn


def make_dirs(d):
    if not os.path.isdir(d):
        log.debug('mkdir -p {!r}'.format(d))
        os.makedirs(d)


def touch(*paths):
    """touch a file.
    """
    msg = 'touch {!r}'.format(paths)
    log.debug(msg)
    for path in paths:
        if os.path.exists(path):
            os.utime(path, None)
        else:
            open(path, 'a').close()
