#!/usr/bin/env python2
import sys, os, re

# A helper script that generates automatic Makefile dependencies from gcc output.

r = re.compile(r"\S+\.(c|h)pp")

def partition(s, separator):
    """Reimplementation of string.partition for Python < 2.5"""
    i = s.find(separator)
    if i < 0:
        return (s,"","")
    return (s[0:i], s[i], s[i+1:])

def sub_func(m):
    filename = m.group(0)
    dir, base = os.path.split(filename)
    base, ext = os.path.splitext(base)
    if ext == ".cpp" or (ext == ".hpp" and os.path.exists("%s/%s.cpp" % (dir, base))):
        return "obj/%s.o" % base

input = sys.stdin.read()                    # Read input from stdin
input = input.replace("\\\n", "")           # Put all in one line
target, sep, depline = partition(input, ":") # Separate target and dependencies
depline = re.sub(r, sub_func, depline)      # Replace dependencies to source files for dependencies on object files
deps = depline.split()                      # Split dependency line into individual dependencies
deps = list(set(deps))                      # Remove duplicates
if target in deps:                          # Remove circular dependencies
    deps.remove(target)
print "%s: %s" % (target, " ".join(deps))   # Join all together and print
