import os
import sys

command = 'gcc Main.c MaxPerm.c -I/usr/local/include/igraph  -L/usr/local/lib -ligraph -lm -o MaxPerm'
os.system(command)

commamd = 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib'
os.system(command)

filename = sys.argv[1]
command = './MaxPerm < ' + filename
os.system(command)

