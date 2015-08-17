VLROOT = ../../common/vlfeat-0.9.18

CXXFLAGS = -std=c++0x -O3 -g -Wall -fmessage-length=0

LIBS_EXTRACTOR = -L$(VLROOT)/bin/glnxa64/ -lvl `pkg-config --cflags --libs opencv`

OBJS_FEATURE_EXTRACTOR = extract.o sift_extractor.o ../../common/feature_set/feature_set.o

OBJS_TEST_FEATURE_EXTRACTOR = test_extract.o sift_extractor.o ../../common/feature_set/feature_set.o

OBJS_TEST_READ = test_read.o ../../common/feature_set/feature_set.o ../../common/file_io/file_io.o

FEATURE_EXTRACTOR = extract

TEST_FEATURE_EXTRACTOR = test_extract

TEST_READ = test_read

all: $(FEATURE_EXTRACTOR) $(TEST_FEATURE_EXTRACTOR) $(TEST_READ)

../../common/file_io/file_io.o: ../../common/file_io/file_io.cc ../../common/file_io/file_io.h
	$(CXX) $(CXXFLAGS) -c ../../common/file_io/file_io.cc \
	-I$(VLROOT) \
	-o ../../common/file_io/file_io.o

../../common/feature_set/feature_set.o: ../../common/feature_set/feature_set.cc ../../common/feature_set/feature_set.h
	$(CXX) $(CXXFLAGS) -c ../../common/feature_set/feature_set.cc \
	-I$(VLROOT) \
	-o ../../common/feature_set/feature_set.o

sift_extractor.o: sift_extractor.cc sift_extractor.h
	$(CXX) $(CXXFLAGS) -c sift_extractor.cc -I$(VLROOT) \
	-o sift_extractor.o $(LIBS_EXTRACTOR)

extract.o: extract.cc sift_extractor.o ../../common/feature_set/feature_set.o
	$(CXX) $(CXXFLAGS) -c extract.cc \
	-I$(VLROOT) \
	-o extract.o $(LIBS_EXTRACTOR)

test_extract.o: test_extract.cc sift_extractor.o ../../common/feature_set/feature_set.o
	$(CXX) $(CXXFLAGS) -c test_extract.cc \
	-I$(VLROOT) \
	-o test_extract.o $(LIBS_EXTRACTOR)

$(FEATURE_EXTRACTOR): $(OBJS_FEATURE_EXTRACTOR)
	$(CXX) $(CXXFLAGS) $(OBJS_FEATURE_EXTRACTOR) \
	-I$(VLROOT) \
	-o $(FEATURE_EXTRACTOR)	$(LIBS_EXTRACTOR)

$(TEST_FEATURE_EXTRACTOR): $(OBJS_TEST_FEATURE_EXTRACTOR)
	$(CXX) $(CXXFLAGS) $(OBJS_TEST_FEATURE_EXTRACTOR) -I$(VLROOT) \
	-o $(TEST_FEATURE_EXTRACTOR) $(LIBS_EXTRACTOR)

$(TEST_READ): $(OBJS_TEST_READ)
	$(CXX) $(CXXFLAGS) $(OBJS_TEST_READ) -I$(VLROOT) \
	-o $(TEST_READ)

clean:
	rm -f $(FEATURE_EXTRACTOR) $(OBJS_FEATURE_EXTRACTOR) \
	$(TEST_FEATURE_EXTRACTOR) $(OBJS_TEST_FEATURE_EXTRACTOR) \
	$(TEST_READ) $(OBJS_TEST_READ) \
	sift_extractor.o
