all: MaxPerm

MaxPerm: Main.c MaxPerm.c
	gcc Main.c MaxPerm.c -I/usr/local/include/igraph  -L/usr/local/lib -ligraph -lm -o MaxPerm

clear:
	rm MaxPerm

network:
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
	./MaxPerm < network10000.dat

email:
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
	./MaxPerm < email-Eu-core.txt

facebook:
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
	./MaxPerm < facebook_combined.txt

wiki:
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
	./MaxPerm < WikiTalk.txt

test:
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
	./MaxPerm < mynetWork.txt

p2p-Gnutella31: MaxPerm
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
	./MaxPerm < p2p-Gnutella31.txt

twitter:
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
	./MaxPerm < twitter.txt

amazon:
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
	./MaxPerm < amazon.txt