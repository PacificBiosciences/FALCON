from __future__ import absolute_import
from __future__ import unicode_literals

from builtins import next
from builtins import range
from builtins import object
from os.path import abspath, expanduser
import io
import contextlib
import gzip
import md5
import re
import subprocess

##
# Utility functions for FastaReader
##


def wrap(s, columns):
    return "\n".join(s[start:start + columns]
                     for start in range(0, len(s), columns))


def splitFastaHeader(name):
    """
    Split a FASTA/FASTQ header into its id and metadata components
    """
    nameParts = re.split('\s', name, maxsplit=1)
    id_ = nameParts[0]
    if len(nameParts) > 1:
        metadata = nameParts[1].strip()
    else:
        metadata = None
    return (id_, metadata)


def splitFileContents(f, delimiter, BLOCKSIZE=8192):
    """
    Same semantics as f.read().split(delimiter), but with memory usage
    determined by largest chunk rather than entire file size
    """
    remainder = io.StringIO()
    while True:
        block = f.read(BLOCKSIZE)
        if not block:
            break
        parts = block.split(delimiter)
        remainder.write(parts[0])
        for part in parts[1:]:
            yield remainder.getvalue()
            remainder = io.StringIO()
            remainder.write(part)
    yield remainder.getvalue()


class FastaRecord(object):
    """
    A FastaRecord object models a named sequence in a FASTA file.
    """
    DELIMITER = ">"
    COLUMNS = 60

    def __init__(self, name, sequence):
        try:
            assert "\n" not in name
            assert "\n" not in sequence
            assert self.DELIMITER not in sequence
            self._name = name
            self._sequence = sequence
            self._md5 = md5.md5(self.sequence).hexdigest()
            self._id, self._metadata = splitFastaHeader(name)
        except AssertionError:
            raise ValueError("Invalid FASTA record data")

    @property
    def name(self):
        """
        The name of the sequence in the FASTA file, equal to the entire
        FASTA header following the '>' character
        """
        return self._name

    @property
    def id(self):
        """
        The id of the sequence in the FASTA file, equal to the FASTA header
        up to the first whitespace.
        """
        return self._id

    @property
    def metadata(self):
        """
        The metadata associated with the sequence in the FASTA file, equal to
        the contents of the FASTA header following the first whitespace
        """
        return self._metadata

    @property
    def sequence(self):
        """
        The sequence for the record as present in the FASTA file.
        (Newlines are removed but otherwise no sequence normalization
        is performed).
        """
        return self._sequence

    @property
    def length(self):
        """
        Get the length of the FASTA sequence
        """
        return len(self._sequence)

    @property
    def md5(self):
        """
        The MD5 checksum (hex digest) of `sequence`
        """
        return self._md5

    @classmethod
    def fromString(cls, s):
        """
        Interprets a string as a FASTA record.  Does not make any
        assumptions about wrapping of the sequence string.
        """
        try:
            lines = s.splitlines()
            assert len(lines) > 1
            assert lines[0][0] == cls.DELIMITER
            name = lines[0][1:]
            sequence = "".join(lines[1:])
            return FastaRecord(name, sequence)
        except AssertionError:
            raise ValueError("String not recognized as a valid FASTA record")

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.name == other.name and
                    self.sequence == other.sequence)
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        """
        Output a string representation of this FASTA record, observing
        standard conventions about sequence wrapping.
        """
        return (">%s\n" % self.name) + \
            wrap(self.sequence, self.COLUMNS)


# These are refactored from ReaderBase/FastaReader.

def yield_fasta_records(f, fn):
    """
    f: fileobj
    fn: str - filename (for exceptions)
    """
    try:
        parts = splitFileContents(f, ">")
        assert "" == next(parts)
        for part in parts:
            yield FastaRecord.fromString(">" + part)
    except AssertionError:
        raise Exception("Invalid FASTA file {!r}".format(fn))


def stream_stdout(call, fn):
    args = call.split()
    proc = subprocess.Popen(args, stdin=open(fn), stdout=subprocess.PIPE)
    return proc.stdout


@contextlib.contextmanager
def open_fasta_reader(fn):
    """
    fn: str - filename

    Note: If you already have a fileobj, you can iterate over yield_fasta_records() directly.

    Streaming reader for FASTA files, useable as a one-shot iterator
    over FastaRecord objects.  Agnostic about line wrapping.
    Example:
    .. doctest::
        TODO: Get data.
        > from pbcore import data
        > filename = data.getTinyFasta()
        > r = FastaReader(filename)
        > with open_fasta_reader(filename) as r:
        ...  for record in r:
        ...     print record.name, len(record.sequence), record.md5
        ref000001|EGFR_Exon_2 183 e3912e9ceacd6538ede8c1b2adda7423
        ref000002|EGFR_Exon_3 203 4bf218da37175a91869033024ac8f9e9
        ref000003|EGFR_Exon_4 215 245bc7a046aad0788c22b071ed210f4d
        ref000004|EGFR_Exon_5 157 c368b8191164a9d6ab76fd328e2803ca
    """
    filename = abspath(expanduser(fn))
    mode = 'r'
    if filename.endswith(".gz"):
        ofs = gzip.open(filename, mode)
    elif filename.endswith(".dexta"):
        ofs = stream_stdout("undexta -vkU -w60 -i", filename)
    else:
        ofs = open(filename, mode)
    yield yield_fasta_records(ofs, filename)
    ofs.close()


class FastaReader(object):
    """Deprecated, but should still work (with filenames).
    """

    def __iter__(self):
        with open_fasta_reader(self.filename) as reader:
            for rec in reader:
                yield rec

    def __init__(self, f):
        self.filename = f
