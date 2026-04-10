# =======================
# Compiler & flags
# =======================
CXX       := g++
CXXFLAGS  := -std=c++17 -O2 -Wall -Wextra
INCLUDES  := -Iinclude -I/usr/include 
LIBS      := -lgmp -lmpfr

# =======================
# Directories
# =======================
SRC_DIR   := src
BUILD_DIR := build

# =======================
# Files
# =======================
MAIN_SRC  := main.cpp
TARGET    := main

# Find all .cpp in src/ recursively
SRC_FILES := $(shell find $(SRC_DIR) -name "*.cpp")

# Object files go to build/, keep folder structure
OBJ_FILES := $(SRC_FILES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
OBJ_FILES += $(BUILD_DIR)/main.o

# =======================
# Rules
# =======================
.PHONY: all clean run

all: $(TARGET)

# Link
$(TARGET): $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

# Compile main.cpp
$(BUILD_DIR)/main.o: $(MAIN_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Compile src/*.cpp → build/*.o
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Run program
run: $(TARGET)
	./$(TARGET)

# Clean
clean:
	rm -rf $(BUILD_DIR) $(TARGET) 
# 	find . -maxdepth 1 -name "*.json" ! -name "trajectory_motion.json" -delete

