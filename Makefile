# compiler
CC = g++
CFLAGS = -Wall -fopenmp -g -O3
INCLUDES =
LDFLAGS = -lgsl -lgslcblas -lm

# directories
SRC_DIRS = src src/Hamiltonian src/matrixElementClasses
OBJ_DIRS = $(patsubst src/%, build/src/%, $(SRC_DIRS))
BIN_DIR = build/bin

# files
SRC_FILES = $(wildcard $(foreach dir, $(SRC_DIRS), $(dir)/*.cpp))
OBJ_FILES = $(patsubst src/%.cpp, build/src/%.o, $(SRC_FILES))

TARGET = $(BIN_DIR)/main

.PHONY: all clean print-debug

all: $(TARGET)

# ensure existing directories
$(OBJ_DIRS) $(BIN_DIR):
	mkdir -p $@

# compile .cpp files
build/src/%.o: src/%.cpp | $(OBJ_DIRS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# main program / linking
$(TARGET): $(OBJ_FILES) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(OBJ_FILES) $(LDFLAGS)

clean:
	rm -rf build/src/Hamiltonian/*
	rm -rf build/src/matrixElementClasses/*
	rm -rf build/src/*.o
	rm -rf build/bin/*

print-debug:
	@echo "SRC_FILES: $(SRC_FILES)"
	@echo "OBJ_FILES: $(OBJ_FILES)"
	@echo "OBJ_DIRS: $(OBJ_DIR)"
	@echo "BIN_DIR: $(BIN_DIR)"

