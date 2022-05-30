CC=g++ --std=c++11 -no-pie

CPPFLAGS = -g -Wall -O3

SLHSD: main.o contig.o aligningFromBam.o scaffoldgraph.o scaffolding.o
	$(CC) -o $@ $^ ./lp/liblpsolve55.a -I $(BAMTOOLS_HOME_INCLUDE)/ $(BAMTOOLS_HOME_LIB)/libbamtools.a -lm -ldl -lz

main.o: main.cpp contig.h aligningFromBam.h scaffoldgraph.h scaffolding.h
	$(CC) -c main.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ -lz
	
contig.o: contig.cpp contig.h
	$(CC) -c contig.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ -lz

aligningFromBam.o: aligningFromBam.cpp contig.h aligningFromBam.h
	$(CC) -c aligningFromBam.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ -lz

scaffoldgraph.o: scaffoldgraph.cpp contig.h scaffoldgraph.h
	$(CC) -c scaffoldgraph.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ -lz

scaffolding.o: scaffolding.cpp contig.h scaffoldgraph.h scaffolding.h
	$(CC) -c scaffolding.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ -lz
	
all: SLHSD
	rm -f *.o
	
clean:
	rm -f *.o
	rm SLHSD
 
 #/home/guanting/anaconda3/envs/bioinfo/include/bamtools/
 #/home/guanting/anaconda3/envs/bioinfo/lib/
 