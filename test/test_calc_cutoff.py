import falcon_kit.mains.calc_cutoff as mod
import helpers
import os.path
import pytest


def test_help():
    try:
        mod.main(['prog', '--help'])
    except SystemExit:
        pass

# Note: genome_size==1 makes math easy.


def test_calc_cutoff(capsys):
    partial_capture_fn = os.path.join(
        helpers.get_test_data_dir(), 'calc_cutoff/partial_capture.txt')
    assert os.path.exists(partial_capture_fn)
    mod.main('prog --coverage 14 1 {}'.format(partial_capture_fn).split())
    out, err = capsys.readouterr()
    assert out == '2'
    assert not err


expected_err = """
GenomeCoverageError: Not enough reads available for desired genome coverage (bases needed=23 > actual=22)
User-provided genome_size: 1
Desired coverage: 23.0
"""


def test_calc_cutoff_err():
    partial_capture_fn = os.path.join(
        helpers.get_test_data_dir(), 'calc_cutoff/partial_capture.txt')
    assert os.path.exists(partial_capture_fn)
    with pytest.raises(Exception) as excinfo:
        mod.main('prog --coverage 23 1 {}'.format(partial_capture_fn).split())
    assert expected_err in str(excinfo.value)


def test_calc_cutoff_errfile(monkeypatch, tmpdir):
    fn = str(tmpdir.mkdir('tmp').join('errfile'))
    monkeypatch.setenv('PBFALCON_ERRFILE', fn)
    partial_capture_fn = os.path.join(
        helpers.get_test_data_dir(), 'calc_cutoff/partial_capture.txt')
    assert os.path.exists(partial_capture_fn)
    with pytest.raises(Exception) as excinfo:
        mod.main('prog --coverage 23 1 {}'.format(partial_capture_fn).split())
    assert expected_err in str(excinfo.value)
    assert expected_err in open(fn).read()
