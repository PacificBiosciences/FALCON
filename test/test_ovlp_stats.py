import falcon_kit.mains.ovlp_stats as mod

def test_help():
    try:
        mod.main(*['prog', '--help'])
    except SystemExit:
        pass
