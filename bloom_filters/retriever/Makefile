CXX=g++
CXXFLAGS = -std=c++0x -O3 -g -Wall -fmessage-length=0

all:

../hasher.o: ../hasher.h ../hasher.cc
	$(CXX) $(CXXFLAGS) -c ../hasher.cc -o ../hasher.o

../bloom.o: ../bloom.h ../bloom.cc
	$(CXX) $(CXXFLAGS) -c ../bloom.cc -o ../bloom.o

../bfindex.o: ../bfindex.h ../bfindex.cc
	$(CXX) $(CXXFLAGS) -c ../bfindex.cc -o ../bfindex.o

clean:
	rm -f ../hasher.o
