#!/usr/bin/env python
import sys, os, h5py, subprocess

def main():
    if len(sys.argv) == 1:
        print "Deletes state data (the wave functions) from itp2d data files, while preserving other data.\nUsage: %s <datafile.h5, datafile2.h5, ...>" % os.path.basename(sys.argv[0])
        sys.exit(1)
    for filename in sys.argv[1:]:
        print "Deleting state data from %s..." % filename
        file = h5py.File(filename)
        del file["/states"]
        file.close()
        outfilename = filename+".stripped"
        subprocess.check_call(["h5repack", "-i", filename, "-o", outfilename])
        print "Done. Stripped datafile saved to %s" % outfilename

if __name__ == "__main__":
    main()
