CXX=icpc

CXXFLAGS = -O3 -qopenmp -std=c++11 -Fast-math -DCPU   -xSSE3 -g -lz
LDFLAGS=-lz

INC_PARAMS=-Ilib -Ilib_cell
srcs = $(wildcard lib_cell/*.cpp)
BUILD_DIR=build

objs := $(srcs:%=$(BUILD_DIR)/%.o)
deps := $(objs:.o=.d)

pace_2: pace_2.cpp $(objs)  
	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@

$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(INC_PARAMS)   $(CXXFLAGS)  -MMD -MP -c $< -o $@ 

Atria_2D_MPI_Ghost_for_0DCell: Atria_2D_MPI_Ghost_for_0DCell.cpp $(objs)  
	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@ $(LDFLAGS)

.PHONY: clean

clean:
	$(RM) $(objs) $(deps) pace

-include $(deps)
MKDIR_P ?= mkdir -p