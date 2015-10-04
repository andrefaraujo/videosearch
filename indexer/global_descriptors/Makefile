YAEL_DIR = ../../common/yael_v260_modif/yael
FEATURESET_DIR = ../../common/feature_set
FILEIO_DIR = ../../common/file_io

CXX=g++
CXXFLAGS = -std=c++0x -O3 -g -Wall -fmessage-length=0 -fopenmp -Wl,-rpath,$(YAEL_DIR)

LIBS_INDEXER = -L$(YAEL_DIR) -lyael

INDEX_DATASET = index_dataset

all: $(INDEX_DATASET)

$(FILEIO_DIR)/file_io.o: $(FILEIO_DIR)/file_io.cc $(FILEIO_DIR)/file_io.h
	$(CXX) $(CXXFLAGS) -c $(FILEIO_DIR)/file_io.cc \
	-o $(FILEIO_DIR)/file_io.o

$(FEATURESET_DIR)/feature_set.o: $(FEATURESET_DIR)/feature_set.cc $(FEATURESET_DIR)/feature_set.h
	$(CXX) $(CXXFLAGS) -c $(FEATURESET_DIR)/feature_set.cc \
	-o $(FEATURESET_DIR)/feature_set.o

gdindex.o: gdindex.h gdindex.cc \
	$(FEATURESET_DIR)/feature_set.o $(FILEIO_DIR)/file_io.o $(YAEL_DIR)/libyael.so
	$(CXX) $(CXXFLAGS) -c gdindex.cc -o gdindex.o \
	-L$(YAEL_DIR) -lyael

$(INDEX_DATASET).o: $(INDEX_DATASET).cc gdindex.o \
	$(YAEL_DIR)/libyael.so
	$(CXX) $(CXXFLAGS) -c $(INDEX_DATASET).cc \
	-I$(YAEL_DIR) \
	-o $(INDEX_DATASET).o -L$(YAEL_DIR) -lyael

$(INDEX_DATASET): $(INDEX_DATASET).o
	$(CXX) $(CXXFLAGS) \
	-I$(YAEL_DIR) \
	-o $(INDEX_DATASET) $(INDEX_DATASET).cc gdindex.o \
	$(FILEIO_DIR)/file_io.o $(FEATURESET_DIR)/feature_set.o $(LIBS_INDEXER)

clean:
	rm -f $(INDEX_DATASET) $(INDEX_DATASET).o gdindex.o
