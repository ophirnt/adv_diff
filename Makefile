CXX        = mpicxx
CXXFLAGS   = -O2 -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include
LDFLAGS    = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc

SRC_DIR    = src
SRCS       = $(SRC_DIR)/main.cpp $(SRC_DIR)/Mesh.hpp
OBJS       = $(SRC_DIR)/main.o
TARGET     = main

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(SRC_DIR)/main.o: $(SRC_DIR)/main.cpp $(SRC_DIR)/Mesh.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(SRC_DIR)/*.o $(TARGET)

.PHONY: all clean
