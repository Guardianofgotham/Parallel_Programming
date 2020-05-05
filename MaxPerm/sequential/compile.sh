gcc Main.c MaxPerm.c -I/usr/local/include/igraph  -L/usr/local/lib -ligraph -lm -o MaxPerm
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
./MaxPerm < network10000.dat
