"""I/O utilities
Not specific to FALCON.
"""
import os
import resource
import shlex
import subprocess as sp
import sys
import traceback

def write_nothing(*args):
    """
    To use,
      LOG = noop
    """

def write_with_pid(*args):
    msg = '[%d]%s\n' %(os.getpid(), ' '.join(args))
    sys.stderr.write(msg)

LOG = write_with_pid

def logstats():
    """This is useful 'atexit'.
    """
    LOG('maxrss:%9d' %(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))

def run_func(args):
    """Wrap multiprocessing.Pool calls.
    Usage:
        pool.imap(run_func, [func, arg0, arg1, ...])
    """
    func = args[0]
    try:
        func_name = func.__name__
    except:
        func_name = repr(func) # but since it must be pickle-able, this should never happen.
    args = args[1:]
    try:
        LOG('starting %s(%s)' %(func_name, ', '.join(repr(a) for a in args)))
        logstats()
        ret = func(*args)
        logstats()
        LOG('finished %s(%s)' %(func_name, ', '.join(repr(a) for a in args)))
        return ret
    except Exception:
        LOG(traceback.format_exc())
        LOG('failed %s(%s)' %(func_name, ', '.join(repr(a) for a in args)))
        return
    except: # KeyboardInterrupt, SystemExit
        LOG('interrupted %s(%s)' %(func_name, ', '.join(repr(a) for a in args)))
        return

def syscall(cmd):
    """Return stdout, fully captured.
    Wait for subproc to finish.
    Raise if empty.
    Raise on non-zero exit-code.
    """
    LOG('$ %s >' %cmd)
    output = sp.check_output(shlex.split(cmd))
    if not output:
        msg = '%r failed to produce any output.' %cmd
        LOG('WARNING: %s' %msg)
    return output

def slurplines(cmd):
    return syscall(cmd).splitlines()

def streamlines(cmd):
    """Stream stdout from cmd.
    Let stderr fall through.
    The returned reader will stop yielding when the subproc exits.
    """
    LOG('$ %s |' %cmd)
    proc = sp.Popen(shlex.split(cmd), stdout=sp.PIPE)
    return proc.stdout
