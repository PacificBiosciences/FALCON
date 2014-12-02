#!/usr/bin/env python

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

import pbcore.io

import sys
"""nucmer -maxmatch all_tigs.fa all_tigs.fa -p all_tigs_self >& /dev/null"""
"""show-coords -o -H -T all_tigs_self.delta | grep CONTAINS | awk '$7>96' | awk '{print $9}' | sort -u > all_tigs_duplicated_ids"""

id_to_remove = set()
with open("all_tigs_duplicated_ids") as f:
    for l in f:
        l = l.strip().split("-")
        major, minor = l[:2]
        id_to_remove.add ( (major, minor) )

f = pbcore.io.FastaReader("all_tigs.fa")
with open("a-tigs_nodup.fa", "w") as f_out:
    for r in f:
        major, minor = r.name.split()[0].split("-")[:2]
        if minor == "0000":
            continue
        if (major, minor) in id_to_remove:
            continue
        if len(r.sequence) < 500:
            continue
        print >>f_out, ">"+r.name
        print >>f_out, r.sequence

f = pbcore.io.FastaReader("primary_tigs_c.fa")
with open("p-tigs_nodup.fa", "w") as f_out:
    for r in f:
        major, minor = r.name.split()[0].split("_")[:2]
        if (major, "0000") in id_to_remove:
            continue
        if len(r.sequence) < 500:
            continue
        print >>f_out, ">"+r.name
        print >>f_out, r.sequence
