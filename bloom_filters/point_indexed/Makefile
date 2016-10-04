YAEL_DIR = ../../common/yael_v260_modif/yael
FEATURESET_DIR = ../../common/feature_set
FILEIO_DIR = ../../common/file_io
GDINDEX_DIR = ../../indexer/global_descriptors

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

all: $(INDEX_DATASET)

$(FILEIO_DIR)/file_io.o: $(FILEIO_DIR)/file_io.cc $(FILEIO_DIR)/file_io.h
	$(CXX) $(CXXFLAGS) -c $(FILEIO_DIR)/file_io.cc \
	-o $(FILEIO_DIR)/file_io.o

$(FEATURESET_DIR)/feature_set.o: $(FEATURESET_DIR)/feature_set.cc $(FEATURESET_DIR)/feature_set.h
	$(CXX) $(CXXFLAGS) -c $(FEATURESET_DIR)/feature_set.cc \
	-o $(FEATURESET_DIR)/feature_set.o

${GDINDEX_DIR}/gdindex.o: ${GDINDEX_DIR}/gdindex.h ${GDINDEX_DIR}/gdindex.cc
	$(CXX) $(CXXFLAGS) -c ${GDINDEX_DIR}/gdindex.cc -o ${GDINDEX_DIR}/gdindex.o

point_index_io.o: point_index_io.h point_index_io.cc
	$(CXX) $(CXXFLAGS) -c point_index_io.cc -o point_index_io.o

$(INDEX_DATASET).o: $(INDEX_DATASET).cc
	$(CXX) $(CXXFLAGS) -c $(INDEX_DATASET).cc \
	-o $(INDEX_DATASET).o

$(INDEX_DATASET): $(INDEX_DATASET).o ${GDINDEX_DIR}/gdindex.o $(FILEIO_DIR)/file_io.o $(FEATURESET_DIR)/feature_set.o point_index_io.o
	$(CXX) $(CXXFLAGS) $(EXEFLAGS) \
	-o $(INDEX_DATASET) $(INDEX_DATASET).o ${GDINDEX_DIR}/gdindex.o point_index_io.o \
	$(FILEIO_DIR)/file_io.o $(FEATURESET_DIR)/feature_set.o $(LIBS_INDEXER)

clean:
	rm -f $(INDEX_DATASET) $(INDEX_DATASET).o $(JOIN_INDEXES) $(JOIN_INDEXES).o ${GDINDEX_DIR}/gdindex.o $(FILEIO_DIR)/file_io.o $(FEATURESET_DIR)/feature_set.o point_index_io.o
