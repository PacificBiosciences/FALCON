from __future__ import unicode_literals
from builtins import str
import falcon_kit.util.system as mod
from falcon_kit.run_support import cd
from falcon_kit.bash import mkdir
import helpers
import pytest
import os
import shutil


def touchtree(*fns):
    for fn in fns:
        mkdir(os.path.dirname(fn))
        with open(fn, 'w'):
            pass


@pytest.yield_fixture
@pytest.fixture(scope="session")
def dirtree(tmpdir_factory):
    tmpdir = tmpdir_factory.mktemp('mytmp')
    with cd(str(tmpdir)):
        touchtree(*dir_picture.strip().split())
        yield
    shutil.rmtree(str(tmpdir))


def picture(fns):
    return '\n'.join(list(fns) + [''])


dir_picture = """\
find_files/file1.txt
find_files/file3.csv
find_files/level_2/file4.txt
find_files/test_file2.txt
find_files/.dotfile
"""

root_dir = 'find_files'


def test_system_find_files(dirtree):
    # Find everything
    result = picture(mod.find_files(root_dir, '*'))
    expected = """\
find_files/.dotfile
find_files/file1.txt
find_files/file3.csv
find_files/test_file2.txt
find_files/level_2/file4.txt
"""
    assert(result == expected)

    # Test normal pattern
    result = picture(mod.find_files(root_dir, '*.txt'))
    expected = """\
find_files/file1.txt
find_files/test_file2.txt
find_files/level_2/file4.txt
"""
    assert(result == expected)

    # Test empty pattern
    result = picture(mod.find_files(root_dir, ''))
    expected = ''
    assert(result == expected)

    # Test non-existent pattern
    result = picture(mod.find_files(root_dir, 'Wubalubadubdub'))
    expected = ''
    assert(result == expected)

    # Test missing dir
    with pytest.raises(Exception) as excinfo:
        list(mod.find_files('Wubalubadubdub', '*'))

    # Test empty dir and empty pattern
    with pytest.raises(Exception) as excinfo:
        list(mod.find_files('', ''))


def test_make_fofn_abs(tmpdir):
    """
    Paths are expanded relative to resolved directory of the input FOFN.
    But otherwise, symlinks do not really need to be resolved.
    """
    with tmpdir.as_cwd():
        expected = os.path.abspath('link') + '\n'
        touchtree('./actual')
        os.symlink('actual', 'link')
        a_fn = './actual.fofn'
        mod.make_dirs('./subdir')
        i_fn = './subdir/i.fofn'
        o_fn = './subdir/o.fofn'
        os.symlink('../actual.fofn', i_fn)
        with open(i_fn, 'w') as stream:
            stream.write('link\n')
        mod.make_fofn_abs(i_fn, o_fn)
        assert open(o_fn).read() == expected
