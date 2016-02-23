#################################################################################$$
# Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted (subject to the limitations in the
# disclaimer below) provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above
#  copyright notice, this list of conditions and the following
#  disclaimer in the documentation and/or other materials provided
#  with the distribution.
#
#  * Neither the name of Pacific Biosciences nor the names of its
#  contributors may be used to endorse or promote products derived
#  from this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#################################################################################$$

from os.path import abspath, expanduser
from cStringIO import StringIO
import md5
import re

def splitFastaHeader( name ):
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
    remainder = StringIO()
    while True:
        block = f.read(BLOCKSIZE)
        if not block:
            break
        parts = block.split(delimiter)
        remainder.write(parts[0])
        for part in parts[1:]:
            yield remainder.getvalue()
            remainder = StringIO()
            remainder.write(part)
    yield remainder.getvalue()

def isFileLikeObject(o):
    return hasattr(o, "read") and hasattr(o, "write")

def getFileHandle(filenameOrFile, mode="r"):
    """
    Given a filename not ending in ".gz", open the file with the
    appropriate mode.
    Given a filename ending in ".gz", return a filehandle to the
    unzipped stream.
    Given a file object, return it unless the mode is incorrect--in
    that case, raise an exception.
    """
    assert mode in ("r", "w")

    if isinstance(filenameOrFile, basestring):
        filename = abspath(expanduser(filenameOrFile))
        if filename.endswith(".gz"):
            return gzip.open(filename, mode)
        else:
            return open(filename, mode)
    elif isFileLikeObject(filenameOrFile):
        return filenameOrFile
    else:
        raise Exception("Invalid type to getFileHandle")


class ReaderBase(object):
    def __init__(self, f):
        """
        Prepare for iteration through the records in the file
        """
        self.file = getFileHandle(f, "r")

    def close(self):
        """
        Close the underlying file
        """
        self.file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


class FastaRecord(object):
    """
    A FastaRecord object models a named sequence in a FASTA file.
    """
    DELIMITER = ">"
    COLUMNS   = 60

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

    def reverseComplement(self, preserveHeader=False):
        """
        Return a new FastaRecord with the reverse-complemented DNA sequence.
        Optionally, supply a name
        """
        rcSequence = sequences.reverseComplement(self.sequence)
        if preserveHeader:
            return FastaRecord(self.name, rcSequence)
        else:
            rcName = '{0} [revcomp]'.format(self.name.strip())
            return FastaRecord(rcName, rcSequence)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.name     == other.name and
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


class FastaReader(ReaderBase):
    """
    Streaming reader for FASTA files, useable as a one-shot iterator
    over FastaRecord objects.  Agnostic about line wrapping.
    Example:
    .. doctest::
        TODO: Get data.
        > from pbcore import data
        > filename = data.getTinyFasta()
        > r = FastaReader(filename)
        > for record in r:
        ...     print record.name, len(record.sequence), record.md5
        ref000001|EGFR_Exon_2 183 e3912e9ceacd6538ede8c1b2adda7423
        ref000002|EGFR_Exon_3 203 4bf218da37175a91869033024ac8f9e9
        ref000003|EGFR_Exon_4 215 245bc7a046aad0788c22b071ed210f4d
        ref000004|EGFR_Exon_5 157 c368b8191164a9d6ab76fd328e2803ca
        >>> r.close()
    """
    DELIMITER = ">"

    def __iter__(self):
        try:
            parts = splitFileContents(self.file, ">")
            assert "" == next(parts)
            for part in parts:
                yield FastaRecord.fromString(">" + part)
        except AssertionError:
            raise ValueError("Invalid FASTA file")

