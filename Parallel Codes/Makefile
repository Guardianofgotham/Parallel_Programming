include $(HCLIB_ROOT)/../modules/system/inc/hclib_system.pre.mak
include $(HCLIB_ROOT)/include/hclib.mak
include $(HCLIB_ROOT)/../modules/system/inc/hclib_system.post.mak

1: P_Main.cpp version1.cpp
	g++ -std=c++14 P_Main.cpp version1.cpp -I/usr/local/include/igraph  -L/usr/local/lib -ligraph -lm $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o parallelExec $(HCLIB_LDLIBS)

2: P_Main.cpp version2.cpp
	g++ -std=c++14 P_Main.cpp version2.cpp -I/usr/local/include/igraph  -L/usr/local/lib -ligraph -lm $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o parallelExec $(HCLIB_LDLIBS)

3: P_Main.cpp version3.cpp
	g++ -std=c++14 P_Main.cpp version3.cpp -I/usr/local/include/igraph  -L/usr/local/lib -ligraph -lm $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o parallelExec $(HCLIB_LDLIBS)

4: P_Main.cpp version4.cpp
	g++ -std=c++14 P_Main.cpp version4.cpp -I/usr/local/include/igraph  -L/usr/local/lib -ligraph -lm $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o parallelExec $(HCLIB_LDLIBS)


test: test.cpp
	g++ -std=c++14 $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) test.cpp -o test $(HCLIB_LDLIBS)


network: parallelExec
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
	./parallelExec < network10000.dat

facebook: parallelExec
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
	./parallelExec < facebook_combined.txt

wiki: parallelExec
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
	./parallelExec < WikiTalk.txt

p2p-Gnutella31: parallelExec
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
	HCLIB_BIND_THREADS=1 HCLIB_STATS=1 HCLIB_WORKERS=6  ./parallelExec < p2p-Gnutella31.txt