# Makefile to build test/bvh_test.cpp which depends on mc_collision/_bvh.hpp

# Compiler
CXX = $(shell uname | grep -q Darwin && echo clang++ || echo g++)
OPTIMIZATION = -O3
CXXFLAGS =  -std=c++17 -flto -lm -Wall -Wextra $(shell uname | grep -q Darwin && echo  -mcpu=apple-m1 -march=armv8-a+simd || echo -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -mavx2 ) $(OPTIMIZATION)

# Directories
SRCDIR = .
INCDIR = ..
OBJDIR = obj
BINDIR = bin
SYSLB2 =
# Source and target files
SRC = $(SRCDIR)/cxx_test.cpp
OBJ = $(OBJDIR)/cxx_test.o
TARGET = $(BINDIR)/cxx_test

# Create necessary directories if they don't exist
$(shell mkdir -p $(OBJDIR))
$(shell mkdir -p $(BINDIR))

# Build target
all: $(TARGET)

# Compile source file into object file
$(OBJ): $(SRC)  $(INCDIR)/cmmcore/vec.h $(INCDIR)/cmmcore/bvh.hpp $(INCDIR)/cmmcore/nurbs_utils.h $(INCDIR)/cmmcore/nurbs.h $(INCDIR)/cmmcore/ccx.h
	$(CXX)  $(CXXFLAGS) -I$(INCDIR) -c $(SRC)  -o $(OBJ) $(SYSLB2)


# Link object file into executable
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS)  $(OBJ) -o $(TARGET) $(SYSLB2)


# Clean object and executable files
clean:
	rm -r $(OBJDIR)
	rm -r $(BINDIR)


.PHONY: all clean