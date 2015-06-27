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

def reprarg(arg):
    if (isinstance(arg, set) or isinstance(arg, list)
            or isinstance(arg, tuple) or isinstance(arg, dict)):
        if len(arg) > 9:
            return '%s(%d elem)' %(type(arg).__name__, len(arg))
    return repr(arg) 
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
        LOG('starting %s(%s)' %(func_name, ', '.join(reprarg(a) for a in args)))
        logstats()
        ret = func(*args)
        logstats()
        LOG('finished %s(%s)' %(func_name, ', '.join(reprarg(a) for a in args)))
        return ret
    except Exception:
        raise Exception(traceback.format_exc())
    except: # KeyboardInterrupt, SystemExit
        LOG('interrupted %s(%s)' %(func_name, ', '.join(reprarg(a) for a in args)))
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
    Note: We do not detect a failure in the underlying process.
    """
    LOG('$ %s |' %cmd)
    proc = sp.Popen(shlex.split(cmd), stdout=sp.PIPE)
    return proc.stdout

class DataReaderContext(object):
    def readlines(self):
        output = self.data.strip()
        for line in output.splitlines():
            yield line
    def __enter__(self):
        pass
    def __exit__(self, *args):
        self.returncode = 0
    def __init__(self, data):
        self.data = data
class ProcessReaderContext(object):
    """Prefer this to slurplines() or streamlines().
    """
    def __enter__(self):
        self.proc = sp.Popen(shlex.split(self.cmd), stdout=sp.PIPE)
    def __exit__(self, etype, evalue, etb):
        if etype is None:
            self.proc.wait()
        else:
            # Exception was raised in "with-block".
            # We cannot wait on proc b/c it might never finish!
            pass
        self.returncode = self.proc.returncode
        if self.returncode:
            msg = "%r <- %r" %(self.returncode, self.cmd)
            raise Exception(msg)
        del self.proc
    def __init__(self, cmd):
        self.cmd = cmd
class CapturedProcessReaderContext(ProcessReaderContext):
    def readlines(self):
        """Usage:

            cmd = 'ls -l'
            reader = CapturedProcessReaderContext(cmd)
            with reader:
                for line in reader.readlines():
                    print line

        Any exception within the 'with-block' is propagated.
        Otherwise, after all lines are read, if 'cmd' failed, Exception is raised.
        """
        output, _ = self.proc.communicate()
        for line in output.splitlines():
            yield line
class StreamedProcessReaderContext(ProcessReaderContext):
    def readlines(self):
        """Usage:

            cmd = 'ls -l'
            reader = StreamedProcessReaderContext(cmd)
            with reader:
                for line in reader.readlines():
                    print line

        Any exception within the 'with-block' is propagated.
        Otherwise, after all lines are read, if 'cmd' failed, Exception is raised.
        """
        for line in self.proc.stdout:
            yield line

def filesize(fn):
    """In bytes.
    Raise if fn does not exist.
    """
    statinfo = os.stat(fn)
    return statinfo.st_size

def validated_fns(fofn):
    """Return list of filenames from fofn.
    Assert none are empty or non-existent.
    """
    fns = open(fofn).read().strip().split("\n")
    for fn in fns:
        assert fn
        assert os.path.isfile(fn)
        assert filesize(fn)
    return fns
