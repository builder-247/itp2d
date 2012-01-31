#!/bin/bash
# This script downloads Google Test, unpacks it in the directory specified by
# the first argument (or to src if no argument is set), and creates a symlink
# named 'gtest' to that directory. The script assumes you have wget and unzip
# installed.

directory=${1-src}
version=gtest-1.6.0
zipfile=${version}.zip

echo Fetching and unpacking Google Test...
wget -nc http://googletest.googlecode.com/files/${zipfile}
unzip -q ${zipfile} -d ${directory}
( cd ${directory} && ln -s -i ${version} gtest )
echo Done
echo It is now safe to remove ${zipfile}
