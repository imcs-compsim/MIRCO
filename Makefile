CC=g++
LDIR=/rzhome/nas/compsim/public/lib/Q1_2015/TPL/trilinos-build/lib
LIBS= $(LDIR)/*.a -lmpi -lmpi_cxx -lblas -lkokkoscore -llapack
IDIR=/rzhome/nas/compsim/public/lib/Q1_2015/TPL/trilinos-build/include
MPI=/usr/lib/openmpi/lib
CFLAGS=-std=c++11 -I$(IDIR) -L$(LDIR) -L$(MPI) -fopenmp

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

bem: main.o
	$(CC) -o bem main.o $(CFLAGS) $(LIBS)

.PHONY:clean

clean:
	rm -f *.o *~ bem
