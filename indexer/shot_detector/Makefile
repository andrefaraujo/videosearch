CXXFLAGS = -std=c++0x -O3 -g -Wall -fmessage-length=0

LIBS = `pkg-config --cflags --libs opencv`

OBJS = main.o

SHOT_DETECTOR = shot_detector

all: $(SHOT_DETECTOR)

main.o: main.cc
	$(CXX) $(CXXFLAGS) -c main.cc \
	-o main.o

$(SHOT_DETECTOR): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) \
	-o $(SHOT_DETECTOR) $(LIBS)

clean:
	rm -f $(SHOT_DETECTOR) main.o
