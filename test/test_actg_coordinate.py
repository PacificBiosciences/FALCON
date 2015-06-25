import falcon_kit.mains.actg_coodinate as mod

def test_help():
    try:
        mod.main(*['prog', '--help'])
    except SystemExit:
        pass
