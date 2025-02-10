# compiler
CC = g++
CFLAGS = -Wall -fopenmp -g -O3
INCLUDES =
LDFLAGS = -lgsl -lgslcblas -lm
SRC_DIR = src
OBJ_DIR = build

# directories
SRC_DIR = src
OBJ_DIR = build/src
BIN_DIR = build/bin
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC_FILES))

TARGET = $(BIN_DIR)/main

# main program / linking
$(TARGET): $(OBJ_FILES) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(OBJ_FILES) $(LDFLAGS)

# compile .cpp files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# create folders if they dont exist
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

clean:
	rm -rf build/src/*.o
	rm -rf build/bin/*
