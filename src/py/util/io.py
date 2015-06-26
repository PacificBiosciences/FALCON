"""I/O utilities
Not specific to FALCON.
"""
import os
import resource
import shlex
import subprocess as sp
import sys

LOG = sys.stderr.write

def write_nothing(*args):
    """
    To use,
      LOG = noop
    """

def logstats():
    """This is useful at run 'atexit'.
    """
    LOG('[%d]maxrss:%9d\n' %(os.getpid(), resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))

def run_func(args):
    func = args[0]
    args = args[1:]
    try:
        LOG('starting %r\n' %func)
        logstats()
        ret = func(*args)
        logstats()
        LOG('finished %r\n' %func)
        return ret
    except KeyboardInterrupt: # and SystemExit?
        return

def syscall(cmd):
    """Return stdout, fully captured.
    Wait for subproc to finish.
    Raise if empty.
    Raise on non-zero exit-code.
    """
    LOG('# %s\n' %cmd)
    output = sp.check_output(shlex.split(cmd))
    if not output:
        msg = '%r failed to produce any output.' %cmd
        LOG('WARNING: %s\n' %msg)
    return output

def slurplines(cmd):
    return syscall(cmd).splitlines()

def streamlines(cmd):
    """Stream stdout from cmd.
    Let stderr fall through.
    The returned reader will stop yielding when the subproc exits.
    """
    LOG('$ %s\n' %cmd)
    proc = sp.Popen(shlex.split(cmd), stdout=sp.PIPE)
    return proc.stdout
