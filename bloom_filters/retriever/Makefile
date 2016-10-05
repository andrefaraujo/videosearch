CXX=g++
CXXFLAGS = -std=c++0x -O3 -g -Wall -fmessage-length=0

all:

../hasher.o: ../hasher.h ../hasher.cc
	$(CXX) $(CXXFLAGS) -c ../hasher.cc -o ../hasher.o

clean:
	rm -f ../hasher.o
