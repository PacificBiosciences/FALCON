import falcon_kit.mains.ovlp_to_graph as mod

def test_help():
    try:
        mod.main(['prog', '--help'])
    except SystemExit:
        pass
