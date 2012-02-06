#!/usr/bin/env python
import sys, os, h5py, subprocess, tempfile, shutil

def main():
    if len(sys.argv) == 1:
        print "Deletes state data (the wave functions) from itp2d data files, while preserving other data.\nUsage: %s <datafile.h5, datafile2.h5, ...>" % os.path.basename(sys.argv[0])
        sys.exit(1)
    filenames = sys.argv[1:]
    for filename in filenames:
        print "Deleting state data from %s..." % filename
        outfilename = "_stripped".join(os.path.splitext(filename))
        # create copy of original file to make sure it is not harmed
        # copy-on-write would be quite useful here
        temphandle, tempname = tempfile.mkstemp(suffix=".h5", prefix=os.path.splitext(filename)[0], dir=os.getcwd())
        shutil.copy(filename, tempname)
        f = h5py.File(tempname)
        del f["/states"]
        f.close()
        subprocess.check_call(["h5repack", "-i", tempname, "-o", outfilename])
        os.unlink(tempname)
        print "Done. Stripped datafile saved to %s" % outfilename

if __name__ == "__main__":
    main()
