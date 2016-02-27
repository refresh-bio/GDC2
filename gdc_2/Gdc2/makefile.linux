all: gdc2
	
CC 	= g++
CFLAGS	= -Wall -O3 -m64 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++11 
CLINK	= -lm -O3 -m64 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++11 -lz

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

gdc2: c1stage.o c2stage.o fasta.o hasher.o main.o p1stage.o qsmodel.o queue.o rangecod.o timer.o
	$(CC) $(CLINK) -o gdc2 c1stage.o c2stage.o fasta.o hasher.o main.o p1stage.o qsmodel.o queue.o rangecod.o timer.o ../libs/libaelf64.a 

clean:
	rm gdc2
	rm *.o
	 