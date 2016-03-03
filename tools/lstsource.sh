#!/bin/sh
#
# Simple script to list all .c, .cpp and .h files in a directory tree.
# Useful for generating the $(project).files for QT Creator, which does not
# handle git branch changes too well.
#

find src test pilot -name "*.cpp" -or -name "*.h" -or -name "*.c" 

