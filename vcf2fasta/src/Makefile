CC=g++

CFLAGS  = -O3 -m64 -std=c++0x -pthread 
LFLAGS  = -O3 -m64 -std=c++0x -lpthread

OMPFLAGS = -fopenmp


all: cut-ref processX VCF2FASTA-d VCF2FASTA-h 

processX: processX.cpp
	$(CC) $(CFLAGS) processX.cpp -o $@

cut-ref: cut-ref.cpp
	$(CC) $(CFLAGS) cut-ref.cpp -o $@

VCF2FASTA-d: VCF2FASTA-d.cpp
	$(CC) $(CFLAGS) VCF2FASTA-d.cpp -o $@

VCF2FASTA-h: VCF2FASTA-h.cpp      
	$(CC) $(CFLAGS) VCF2FASTA-h.cpp -o $@  

.cpp.o:
	$(CC) $(CFLAGS) -c $<
clean:
	-rm cut-ref
	-rm processX
	-rm VCF2FASTA-d
	-rm VCF2FASTA-h
