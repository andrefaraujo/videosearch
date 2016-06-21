YAEL_DIR = ../common/yael_v260_modif/yael
GDINDEX_DIR = ../indexer/global_descriptors
FEATURESET_DIR = ../common/feature_set
FILEIO_DIR = ../common/file_io

CXX=g++
CXXFLAGS = -std=c++0x -O3 -g -Wall -fmessage-length=0
EXEFLAGS = -Wl,-rpath,'$$ORIGIN'/$(YAEL_DIR)

LIBS_RETRIEVER = -L$(YAEL_DIR) -lyael

PLATFORM=$(shell uname)
ifeq ($(PLATFORM),Darwin) # Mac options
else # Linux options
	CXXFLAGS += -fopenmp
endif

RETRIEVE_DATASET = retrieve_on_dataset

all: $(RETRIEVE_DATASET)

$(FILEIO_DIR)/file_io.o: $(FILEIO_DIR)/file_io.cc $(FILEIO_DIR)/file_io.h
	$(CXX) $(CXXFLAGS) -c $(FILEIO_DIR)/file_io.cc \
	-o $(FILEIO_DIR)/file_io.o

$(FEATURESET_DIR)/feature_set.o: $(FEATURESET_DIR)/feature_set.cc $(FEATURESET_DIR)/feature_set.h
	$(CXX) $(CXXFLAGS) -c $(FEATURESET_DIR)/feature_set.cc \
	-o $(FEATURESET_DIR)/feature_set.o

$(GDINDEX_DIR)/gdindex.o: $(GDINDEX_DIR)/gdindex.h $(GDINDEX_DIR)/gdindex.cc
	$(CXX) $(CXXFLAGS) -c $(GDINDEX_DIR)/gdindex.cc -o $(GDINDEX_DIR)/gdindex.o

retriever.o: retriever.cc retriever.h
	$(CXX) $(CXXFLAGS) -c retriever.cc \
	-o retriever.o

$(RETRIEVE_DATASET).o: $(RETRIEVE_DATASET).cc
	$(CXX) $(CXXFLAGS) -c $(RETRIEVE_DATASET).cc \
	-o $(RETRIEVE_DATASET).o

$(RETRIEVE_DATASET): $(RETRIEVE_DATASET).o retriever.o $(FILEIO_DIR)/file_io.o $(FEATURESET_DIR)/feature_set.o $(GDINDEX_DIR)/gdindex.o
	$(CXX) $(CXXFLAGS) $(EXEFLAGS) \
	-o $(RETRIEVE_DATASET) $(RETRIEVE_DATASET).o retriever.o $(GDINDEX_DIR)/gdindex.o \
	$(FILEIO_DIR)/file_io.o $(FEATURESET_DIR)/feature_set.o $(LIBS_RETRIEVER)

clean:
	rm -f $(RETRIEVE_DATASET) $(RETRIEVE_DATASET).o retriever.o $(FILEIO_DIR)/file_io.o $(FEATURESET_DIR)/feature_set.o $(GDINDEX_DIR)/gdindex.o
