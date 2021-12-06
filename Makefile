CC=g++
LDIR=/imcs/public/compsim/lib/Q4_2019_ubuntu20/TPL/trilinos-build-release/lib
LIBS= $(LDIR)/*.a -lmpi -lmpi_cxx -lblas -lkokkoscore -llapack -ldl
IDIR=/imcs/public/compsim/lib/Q4_2019_ubuntu20/TPL/trilinos-build-release/include
MPI=/usr/lib/x86_64-linux-gnu/openmpi/lib
CFLAGS=-std=c++11 -I$(IDIR) -L$(LDIR) -L$(MPI) -fopenmp -ljsoncpp

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

bem: main.o topology.o
	$(CC) -o bem main.o topology.o $(CFLAGS) $(LIBS)

.PHONY:clean

clean:
	rm -f *.o *~ bem
