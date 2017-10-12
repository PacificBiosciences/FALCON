import falcon_kit.mains.ovlp_to_graph as mod
import pytest


def test_help():
    try:
        mod.main(['prog', '--help'])
    except SystemExit:
        pass


def test_reverse_end():
    ret = mod.reverse_end('123:B')
    expected = '123:E'
    assert(expected == ret)

    ret = mod.reverse_end('456:E')
    expected = '456:B'
    assert(expected == ret)

    ret = mod.reverse_end(':B')
    expected = ':E'
    assert(expected == ret)

    ret = mod.reverse_end(':E')
    expected = ':B'
    assert(expected == ret)

    with pytest.raises(Exception) as e_info:
        ret = mod.reverse_end(':')

    with pytest.raises(Exception) as e_info:
        ret = mod.reverse_end('')

    with pytest.raises(Exception) as e_info:
        ret = mod.reverse_end('12345:C')

    with pytest.raises(Exception) as e_info:
        ret = mod.reverse_end(':::')
