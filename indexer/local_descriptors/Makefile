VLROOT = ../../common/vlfeat-0.9.18
FEATURESET_DIR = ../../common/feature_set
FILEIO_DIR = ../../common/file_io

FEATURE_EXTRACTOR = extract

TEST_FEATURE_EXTRACTOR = test_extract

TEST_READ = test_read

# Define mac commands for linking
MAC_INSTALL_CMD=install_name_tool -change "@loader_path/libvl.dylib" "$(VLROOT)/bin/maci64/libvl.dylib"
PLATFORM=$(shell uname)
ifeq ($(PLATFORM),Darwin)
	VLROOT_ARCH = maci64
	MAC_INSTALL_CMD_FE = $(MAC_INSTALL_CMD) $(FEATURE_EXTRACTOR)
	MAC_INSTALL_CMD_TFE = $(MAC_INSTALL_CMD) $(TEST_FEATURE_EXTRACTOR)
	MAC_INSTALL_CMD_TR = $(MAC_INSTALL_CMD) $(TEST_READ)
else
	VLROOT_ARCH = glnxa64
endif

CXXFLAGS = -std=c++0x -O3 -g -Wall -fmessage-length=0 -Wl,-rpath,$(VLROOT)/bin/$(VLROOT_ARCH)

LIBS_EXTRACTOR = -L$(VLROOT)/bin/$(VLROOT_ARCH)/ -lvl `pkg-config --cflags --libs opencv`

OBJS_FEATURE_EXTRACTOR = sift_extractor.o $(FEATURESET_DIR)/feature_set.o

OBJS_TEST_FEATURE_EXTRACTOR = sift_extractor.o $(FEATURESET_DIR)/feature_set.o

OBJS_TEST_READ = $(FEATURESET_DIR)/feature_set.o $(FILEIO_DIR)/file_io.o


all: $(FEATURE_EXTRACTOR) $(TEST_FEATURE_EXTRACTOR) $(TEST_READ)

$(FILEIO_DIR)/file_io.o: $(FILEIO_DIR)/file_io.cc $(FILEIO_DIR)/file_io.h
	$(CXX) $(CXXFLAGS) -c $(FILEIO_DIR)/file_io.cc \
	-o $(FILEIO_DIR)/file_io.o

$(FEATURESET_DIR)/feature_set.o: $(FEATURESET_DIR)/feature_set.cc $(FEATURESET_DIR)/feature_set.h
	$(CXX) $(CXXFLAGS) -c $(FEATURESET_DIR)/feature_set.cc \
	-o $(FEATURESET_DIR)/feature_set.o

sift_extractor.o: sift_extractor.cc sift_extractor.h
	$(CXX) $(CXXFLAGS) -c sift_extractor.cc -I$(VLROOT) \
	-o sift_extractor.o $(LIBS_EXTRACTOR)

extract.o: extract.cc $(OBJS_FEATURE_EXTRACTOR)
	$(CXX) $(CXXFLAGS) -c extract.cc \
	-I$(VLROOT) \
	-o extract.o $(LIBS_EXTRACTOR)

test_extract.o: test_extract.cc $(OBJS_TEST_FEATURE_EXTRACTOR)
	$(CXX) $(CXXFLAGS) -c test_extract.cc \
	-I$(VLROOT) \
	-o test_extract.o $(LIBS_EXTRACTOR)

test_read.o: test_read.cc $(OBJS_TEST_READ)
	$(CXX) $(CXXFLAGS) -c test_read.cc \
	-I$(VLROOT) \
	-o test_read.o

$(FEATURE_EXTRACTOR): $(FEATURE_EXTRACTOR).o
	$(CXX) $(CXXFLAGS) $(FEATURE_EXTRACTOR).o $(OBJS_FEATURE_EXTRACTOR) \
	-I$(VLROOT) \
	-o $(FEATURE_EXTRACTOR)	$(LIBS_EXTRACTOR); \
	$(MAC_INSTALL_CMD_FE)

$(TEST_FEATURE_EXTRACTOR): $(TEST_FEATURE_EXTRACTOR).o
	$(CXX) $(CXXFLAGS) $(TEST_FEATURE_EXTRACTOR).o $(OBJS_TEST_FEATURE_EXTRACTOR) -I$(VLROOT) \
	-o $(TEST_FEATURE_EXTRACTOR) $(LIBS_EXTRACTOR); \
	$(MAC_INSTALL_CMD_TFE)

$(TEST_READ): $(TEST_READ).o
	$(CXX) $(CXXFLAGS) $(TEST_READ).o $(OBJS_TEST_READ) -I$(VLROOT) \
	-o $(TEST_READ); \
	$(MAC_INSTALL_CMD_TR)

clean:
	rm -f $(FEATURE_EXTRACTOR) $(FEATURE_EXTRACTOR).o $(OBJS_FEATURE_EXTRACTOR) \
	$(TEST_FEATURE_EXTRACTOR) $(TEST_FEATURE_EXTRACTOR).o $(OBJS_TEST_FEATURE_EXTRACTOR) \
	$(TEST_READ) $(TEST_READ).o $(OBJS_TEST_READ) \
	sift_extractor.o
