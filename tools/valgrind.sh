#!/bin/sh
vg=valgrind
opts=

case `hostname` in 
	stallo*)
	export MPIWRAP_DEBUG=verbose
	export \
	LD_PRELOAD=/global/apps/valgrind/lib/valgrind/libmpiwrap-amd64-linux.so 
	;;
esac

while true; do
	[ $# = 0 ] && break
	if [ $1 = "leak" ]; then
		opts="$opts --leak-check=full"
	elif [ $1 = "reach" ]; then
		opts="$opts --leak-check=full --show-reachable=yes"
	elif [ $1 = "track" ]; then
		opts="$opts --track-origins=yes" 
	elif [ $1 = "omp" ]; then
		opts="$opts --tool=helgrind" 
	elif [ $1 = "mpi" ]; then
		vg="mpirun -np 2 $vg"
	fi
	shift
done

[ `uname` = "Darwin" ] && dsymutil src/mrchem.bin 

echo "$vg $opts src/mrchem.bin 2>&1 |tee valgrind.out"
$vg $opts src/mrchem.bin 2>&1 |tee valgrind.out

