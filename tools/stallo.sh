#!/bin/sh

topdir=$HOME/src/mrchem

cd $topdir
rsync -avc \
--exclude '.*' \
--exclude 'share/*' \
--exclude 'external/*' \
--exclude 'Documentation/*' \
--exclude 'Debug/*' \
--exclude 'Release/*' \
. stallo:src/mrchem

#cd $topdir
#if [ $# = 1 ]; then
#    ./config/build_remote.sh stallo src init
#else
#    ./config/build_remote.sh stallo src
#fi

