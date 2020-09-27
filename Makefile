CC=g++
LDIR=/imcs/Q1_2015/trilinos/lib
LIBS= $(LDIR)/*.a -lmpi /imcs/lib/libopenblas.a -lkokkoscore /usr/lib64/liblapack.so.3
IDIR=/imcs/Q1_2015/trilinos/include
MPI=/opt/openmpi/4.0.0/gcc
CFLAGS=-std=c++11 -O3 -I$(IDIR) -I$(MPI)/include -L$(LDIR) -L$(MPI)/lib -fopenmp

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

bem: main.o
	$(CC) -o bem main.o -Wl,--copy-dt-needed-entries $(CFLAGS) $(LIBS)

.PHONY:clean

clean:
	rm -f *.o *~ bem
