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

CXXFLAGS = -std=c++0x -O3 -g -Wall -fmessage-length=0
EXEFLAGS = -Wl,-rpath,'$$ORIGIN'/$(VLROOT)/bin/$(VLROOT_ARCH)

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
	$(CXX) $(CXXFLAGS) -c sift_extractor.cc \
	-o sift_extractor.o

$(FEATURE_EXTRACTOR).o: extract.cc
	$(CXX) $(CXXFLAGS) -c extract.cc \
	-o extract.o

$(TEST_FEATURE_EXTRACTOR).o: test_extract.cc
	$(CXX) $(CXXFLAGS) -c test_extract.cc \
	-o test_extract.o

$(TEST_READ).o: test_read.cc
	$(CXX) $(CXXFLAGS) -c test_read.cc \
	-o test_read.o

$(FEATURE_EXTRACTOR): $(FEATURE_EXTRACTOR).o $(OBJS_FEATURE_EXTRACTOR)
	$(CXX) $(CXXFLAGS) $(EXEFLAGS) $(FEATURE_EXTRACTOR).o $(OBJS_FEATURE_EXTRACTOR) \
	-o $(FEATURE_EXTRACTOR)	$(LIBS_EXTRACTOR); \
	$(MAC_INSTALL_CMD_FE)

$(TEST_FEATURE_EXTRACTOR): $(TEST_FEATURE_EXTRACTOR).o $(OBJS_TEST_FEATURE_EXTRACTOR)
	$(CXX) $(CXXFLAGS) $(EXEFLAGS) $(TEST_FEATURE_EXTRACTOR).o $(OBJS_TEST_FEATURE_EXTRACTOR) \
	-o $(TEST_FEATURE_EXTRACTOR) $(LIBS_EXTRACTOR); \
	$(MAC_INSTALL_CMD_TFE)

$(TEST_READ): $(TEST_READ).o $(OBJS_TEST_READ)
	$(CXX) $(CXXFLAGS) $(EXEFLAGS) $(TEST_READ).o $(OBJS_TEST_READ) \
	-o $(TEST_READ); \
	$(MAC_INSTALL_CMD_TR)

clean:
	rm -f $(FEATURE_EXTRACTOR) $(FEATURE_EXTRACTOR).o $(OBJS_FEATURE_EXTRACTOR) \
	$(TEST_FEATURE_EXTRACTOR) $(TEST_FEATURE_EXTRACTOR).o $(OBJS_TEST_FEATURE_EXTRACTOR) \
	$(TEST_READ) $(TEST_READ).o $(OBJS_TEST_READ) \
	sift_extractor.o $(FEATURESET_DIR)/feature_set.o $(FILEIO_DIR)/file_io.o
