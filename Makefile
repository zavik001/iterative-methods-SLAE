CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -pedantic

SRC_DIR = src
INCLUDE_DIR = include
BIN_DIR = bin
DATA_DIR = data

TARGET = $(BIN_DIR)/main

SRC_FILES = $(SRC_DIR)/SLAE.cpp Main.cpp

all: $(TARGET)

$(TARGET): $(SRC_FILES)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) $(SRC_FILES) -o $(TARGET)

run: all
	./$(TARGET)

clean:
	rm -rf $(BIN_DIR)

.PHONY: all run clean
