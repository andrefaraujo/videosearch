YAEL_DIR = ../common/yael_v438_modif/yael
GDINDEX_DIR = ../indexer/global_descriptors
FEATURESET_DIR = ../common/feature_set
FILEIO_DIR = ../common/file_io

CXX=g++
CXXFLAGS = -std=c++0x -O3 -g -Wall -fmessage-length=0 -fopenmp -Wl,-rpath,$(YAEL_DIR)

LIBS_RETRIEVER = -L$(YAEL_DIR) -lyael

RETRIEVE_DATASET = retrieve_on_dataset

all: $(RETRIEVE_DATASET)

$(FILEIO_DIR)/file_io.o: $(FILEIO_DIR)/file_io.cc $(FILEIO_DIR)/file_io.h
	$(CXX) $(CXXFLAGS) -c $(FILEIO_DIR)/file_io.cc \
	-o $(FILEIO_DIR)/file_io.o

$(FEATURESET_DIR)/feature_set.o: $(FEATURESET_DIR)/feature_set.cc $(FEATURESET_DIR)/feature_set.h
	$(CXX) $(CXXFLAGS) -c $(FEATURESET_DIR)/feature_set.cc \
	-o $(FEATURESET_DIR)/feature_set.o

$(GDINDEX_DIR)/gdindex.o: $(GDINDEX_DIR)/gdindex.h $(GDINDEX_DIR)/gdindex.cc \
	$(FEATURESET_DIR)/feature_set.o $(FILEIO_DIR)/file_io.o $(YAEL_DIR)/libyael.so
	$(CXX) $(CXXFLAGS) -c $(GDINDEX_DIR)/gdindex.cc -o $(GDINDEX_DIR)/gdindex.o \
	-L$(YAEL_DIR) -lyael

retriever.o: $(GDINDEX_DIR)/gdindex.o retriever.cc retriever.h
	$(CXX) $(CXXFLAGS) -c retriever.cc \
	-I$(YAEL_DIR) \
	-o retriever.o -L$(YAEL_DIR) -lyael

$(RETRIEVE_DATASET).o: $(RETRIEVE_DATASET).cc retriever.o $(GDINDEX_DIR)/gdindex.o \
	$(YAEL_DIR)/libyael.so
	$(CXX) $(CXXFLAGS) -c $(RETRIEVE_DATASET).cc \
	-I$(YAEL_DIR) \
	-o $(RETRIEVE_DATASET).o -L$(YAEL_DIR) -lyael

$(RETRIEVE_DATASET): $(RETRIEVE_DATASET).o
	$(CXX) $(CXXFLAGS) \
	-I$(YAEL_DIR) \
	-o $(RETRIEVE_DATASET) $(RETRIEVE_DATASET).cc retriever.o $(GDINDEX_DIR)/gdindex.o $(LIBS_RETRIEVER)

clean:
	rm -f $(RETRIEVE_DATASET) $(RETRIEVE_DATASET).o retriever.o
