# Makefile to build test/bvh_test.cpp which depends on mc_collision/_bvh.hpp

# Compiler
CXX = emcc
OPTIMIZATION = -O3
CXXFLAGS = $(OPTIMIZATION) -std=c++17 -s MODULARIZE=1 -s EXPORT_NAME="cmmcore" --bind

# Directories
SRCDIR = .
INCDIR = ..
OBJDIR = .
BINDIR = .
SYSLB2 =
# Source and target files
SRC = $(SRCDIR)/cmmcore.cpp
OBJ = $(OBJDIR)/cmmcore.wasm
TARGET = $(BINDIR)/cmmcore.js

# Create necessary directories if they don't exist
#$(shell mkdir -p $(OBJDIR))
#$(shell mkdir -p $(BINDIR))

# Build target
all: $(TARGET)

# Compile source file into object file
#$(OBJ): $(SRC)
#	$(CXX) -o  $(CXXFLAGS) -I$(INCDIR) -c $(SRC)  -o $(OBJ) $(SYSLB2)


# Compile source file
$(TARGET):
	$(CXX) $(SRC) -o $(TARGET) -I$(INCDIR) $(CXXFLAGS)



# Clean object and executable files
clean:
	rm -rf $(OBJ)
	rm -rf $(TARGET)


.PHONY: all clean