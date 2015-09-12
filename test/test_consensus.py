import falcon_kit.mains.consensus as mod

def test_help():
    try:
        mod.main(['prog', '--help'])
    except SystemExit:
        pass
