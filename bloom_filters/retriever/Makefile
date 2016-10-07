CXX=g++
CXXFLAGS = -std=c++0x -O3 -g -Wall -fmessage-length=0

RETRIEVE_DATASET = retrieve_on_dataset_bf

all: $(RETRIEVE_DATASET)

../hasher.o: ../hasher.h ../hasher.cc
	$(CXX) $(CXXFLAGS) -c ../hasher.cc -o ../hasher.o

../bloom.o: ../bloom.h ../bloom.cc
	$(CXX) $(CXXFLAGS) -c ../bloom.cc -o ../bloom.o

../bfindex.o: ../bfindex.h ../bfindex.cc
	$(CXX) $(CXXFLAGS) -c ../bfindex.cc -o ../bfindex.o

../point_indexed/point_index_io.o: ../point_indexed/point_index_io.h ../point_indexed/point_index_io.cc
	$(CXX) $(CXXFLAGS) -c ../point_indexed/point_index_io.cc -o ../point_indexed/point_index_io.o

$(RETRIEVE_DATASET).o: $(RETRIEVE_DATASET).cc
	$(CXX) $(CXXFLAGS) -c $(RETRIEVE_DATASET).cc \
	-o $(RETRIEVE_DATASET).o

$(RETRIEVE_DATASET): $(RETRIEVE_DATASET).o ../hasher.o ../bloom.o ../bfindex.o ../point_indexed/point_index_io.o
	$(CXX) $(CXXFLAGS) $(EXEFLAGS) \
	-o $(RETRIEVE_DATASET) $(RETRIEVE_DATASET).o ../hasher.o ../bloom.o \
	../bfindex.o ../point_indexed/point_index_io.o

clean:
	rm -f ../hasher.o ../bloom.o ../bfindex.o $(RETRIEVE_DATASET).o $(RETRIEVE_DATASET) ../point_indexed/point_index_io.o
