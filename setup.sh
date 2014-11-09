module add gcc/4.8.1
module add mpich2/1.4

export LIBRARY_PATH=$LIBRARY_PATH:./lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./lib

make

