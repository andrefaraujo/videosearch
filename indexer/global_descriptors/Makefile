YAEL_DIR = ../../common/yael_v260_modif/yael
FEATURESET_DIR = ../../common/feature_set
FILEIO_DIR = ../../common/file_io

CXX=g++
CXXFLAGS = -std=c++0x -O3 -g -Wall -fmessage-length=0
EXEFLAGS = -Wl,-rpath,'$$ORIGIN'/$(YAEL_DIR)

LIBS_INDEXER = -L$(YAEL_DIR) -lyael

PLATFORM=$(shell uname)
ifeq ($(PLATFORM),Darwin) # Mac options
else # Linux options
	CXXFLAGS += -fopenmp
endif

INDEX_DATASET = index_dataset
JOIN_INDEXES = join_indexes

all: $(INDEX_DATASET) $(JOIN_INDEXES)

$(FILEIO_DIR)/file_io.o: $(FILEIO_DIR)/file_io.cc $(FILEIO_DIR)/file_io.h
	$(CXX) $(CXXFLAGS) -c $(FILEIO_DIR)/file_io.cc \
	-o $(FILEIO_DIR)/file_io.o

$(FEATURESET_DIR)/feature_set.o: $(FEATURESET_DIR)/feature_set.cc $(FEATURESET_DIR)/feature_set.h
	$(CXX) $(CXXFLAGS) -c $(FEATURESET_DIR)/feature_set.cc \
	-o $(FEATURESET_DIR)/feature_set.o

gdindex.o: gdindex.h gdindex.cc
	$(CXX) $(CXXFLAGS) -c gdindex.cc -o gdindex.o

$(INDEX_DATASET).o: $(INDEX_DATASET).cc
	$(CXX) $(CXXFLAGS) -c $(INDEX_DATASET).cc \
	-o $(INDEX_DATASET).o

$(INDEX_DATASET): $(INDEX_DATASET).o gdindex.o $(FILEIO_DIR)/file_io.o $(FEATURESET_DIR)/feature_set.o
	$(CXX) $(CXXFLAGS) $(EXEFLAGS) \
	-o $(INDEX_DATASET) $(INDEX_DATASET).o gdindex.o \
	$(FILEIO_DIR)/file_io.o $(FEATURESET_DIR)/feature_set.o $(LIBS_INDEXER)

$(JOIN_INDEXES).o: $(JOIN_INDEXES).cc
	$(CXX) $(CXXFLAGS) -c $(JOIN_INDEXES).cc \
	-o $(JOIN_INDEXES).o

$(JOIN_INDEXES): $(JOIN_INDEXES).o gdindex.o $(FILEIO_DIR)/file_io.o $(FEATURESET_DIR)/feature_set.o
	$(CXX) $(CXXFLAGS) $(EXEFLAGS) \
	-o $(JOIN_INDEXES) $(JOIN_INDEXES).o gdindex.o \
	$(FILEIO_DIR)/file_io.o $(FEATURESET_DIR)/feature_set.o $(LIBS_INDEXER)

clean:
	rm -f $(INDEX_DATASET) $(INDEX_DATASET).o $(JOIN_INDEXES) $(JOIN_INDEXES).o gdindex.o $(FILEIO_DIR)/file_io.o $(FEATURESET_DIR)/feature_set.o 
