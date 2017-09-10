from .io import system
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
    for d, paths in dir2paths.iteritems():
        bases = [os.path.basename(path) for path in paths]
        base2rel = {os.path.basename(path): os.path.relpath(path, d) for path in paths}
        assert len(base2rel) == len(bases), 'Non-unique basename in {}'.format(repr(paths))
        for existing_base in os.listdir(d):
            existing_path = os.path.join(d, existing_base)
            if os.path.islink(existing_path):
                if existing_base in base2rel:
                    if os.readlink(existing_path) != base2rel[existing_base]:
                        os.unlink(existing_path) # Wrong target (or non-relative) so remove it.
                    else:
                        del base2rel[existing_base] # Just keep it.
                else:
                    os.unlink(existing_path) # Old? Remove it for safety.
        for base, rel in base2rel.iteritems():
            path = os.path.join(d, base)
            os.symlink(rel, path)

def lfs_setstripe_maybe(path='.', stripe=12):
    path = os.path.abspath(path)
    rc = system('lfs setstripe -c {:d} {!s}'.format(stripe, path))
    if rc:
        log.info('Apparently {!r} is not lustre in filesystem.'.format(path))
    else:
        log.info('This lfs stripe ({}) should propagate to subdirs of {!r}.'.format(stripe, path))

def find_files(root_path, pattern):
    """
    Finds all files with filenames formatted as the
    given pattern, descending down from root_path.
    """
    for root, dirs, files in os.walk(root_path):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(root, filename)
