CC=g++
LDIR=/imcs/Q1_2015/trilinos/lib
LIBS= $(LDIR)/*.a /imcs/lib/libopenblas.a -lkokkoscore
IDIR=/imcs/Q1_2015/trilinos/include
MPI=/opt/openmpi/4.0.0/gcc/lib
CFLAGS=-std=c++11 -I$(IDIR) -L$(LDIR) -L$(MPI) -fopenmp

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

bem: main.o
	$(CC) -o test_program main.o $(CFLAGS) $(LIBS)

.PHONY:clean

clean:
	rm -f *.o *~ test_program
