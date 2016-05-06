#!/usr/bin/env python

import subprocess
import os.path
import argparse
import time
import fcntl
import sys
import os
import datetime

# DG, May 3, 2016

sys.stderr.write( "from copy_big_files.py\n" )

# alert!  This script is running in the 0-rawreads/preads subdirectory
# of the network disk

ap = argparse.ArgumentParser(description="copy_big_files.py /var/tmp")
ap.add_argument( "szDevShmRawReadsDir" )
args = ap.parse_args()

print "from copy_big_files.py"
sys.stderr.write( "from copy_big_files.py\n" )
szCommand = "uname -a"
szOutput = subprocess.check_output( szCommand, shell = True )
print "from: " + szOutput,
sys.stderr.write( "from: " + szOutput + "\n" )
# looks like:
# > uname -a
#Linux e61.grid.gs.washington.edu 2.6.32-573.18.1.el6.x86_64 #1 SMP Wed Jan 6 11:20:49 EST 2016 x86_64 x86_64 x86_64 GNU/Linux

aWords = szOutput.split()
szNodeFull = aWords[1]
aWords2 = szNodeFull.split(".")
szNode = aWords2[0]


szCommand = "mkdir -p ../locks"
subprocess.call(szCommand, shell = True )

szLockFile = "../locks/copy_big_files_" + szNode + "_lock"

with open( szLockFile, "w" ) as fLockFile:
    # wait until have lock
    print "about to try to get lock"
    date1 = datetime.datetime.now()
    fcntl.flock( fLockFile.fileno(), fcntl.LOCK_EX )
    date2 = datetime.datetime.now()
    
    nSecondsElapsed = ( date2 - date1 ).total_seconds()
    print "{:f} seconds to get lock".format( nSecondsElapsed )
    

    if ( os.path.exists( args.szDevShmRawReadsDir + "/raw_reads.db" ) ):
        print "file " + args.szDevShmRawReadsDir + "/raw_reads.db" + " already exists"
        fcntl.flock( fLockFile.fileno(), fcntl.LOCK_UN )
        sys.exit( 0 )

    szCommand = "mkdir -p " + args.szDevShmRawReadsDir
    subprocess.call( szCommand, shell = True )

    # if reached here, we have the lock and the file doesn't exist

    szBwlimit = "50000"

    aFilesToCopy = [".raw_reads.bps", ".raw_reads.idx", "raw_reads.db" ]

    

    for szFileToCopy in aFilesToCopy:
        date3 = datetime.datetime.now()
        szCommand = "rsync -v --bwlimit=" + szBwlimit + " ../" + szFileToCopy + " " + args.szDevShmRawReadsDir

        sys.stderr.write( "about to " + szCommand + "\n" )
        subprocess.call( szCommand, shell = True )
        date4 = datetime.datetime.now()
        fSecondsElapsed = ( date4 - date3 ).total_seconds()

        ( nMinutes, fSeconds ) = divmod( nSecondsElapsed, 60 )

        
        print "completed copying {:s} in {:f}:{:02.0f}".format( szFileToCopy, nMinutes, fSeconds )

        szOutput =  "completed copying {:s} in {:f}:{:02.0f}".format( szFileToCopy, nMinutes, fSeconds )
        sys.stderr.write( szOutput + "\n" )


    fcntl.flock( fLockFile.fileno(), fcntl.LOCK_UN )
    sys.exit( 0 )













