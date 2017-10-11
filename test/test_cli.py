import commands


def test_gen_gfa_v1():
    cmd = 'python -m falcon_kit.mains.gen_gfa_v1 --help'
    rc, out = commands.getstatusoutput(cmd)
    assert rc == 0
    assert out != ''
