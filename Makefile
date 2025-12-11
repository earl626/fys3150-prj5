
SRC_DIR = src
BUILD_DIR = build
OUTPUT_DIR = output
DIRS = $(BUILD_DIR) $(OUTPUT_DIR) # List of directories that must exist

CXX = g++
FLAGS = -std=c++17 -O2 -Wall -Wextra -Wpedantic
FLAGS_OPEN_MP = -std=c++17 -O3 -fopenmp -Wall -Wextra -Wpedantic
INCLUDE = -I include -larmadillo

# Create directories if they do not exist
$(DIRS):
	mkdir -p $@

SRC_MAIN = $(SRC_DIR)/main.cpp $(SRC_DIR)/utils.cpp $(SRC_DIR)/solver.cpp
main: $(DIRS) $(SRC_MAIN)
	$(CXX) $(FLAGS) $(SRC_MAIN) -o $(BUILD_DIR)/main $(INCLUDE)
