#CXXFLAGS = -O3 -std=c++17 -Iinclude -Xpreprocessor -fopenmp -lomp -I"$(shell brew --prefix libomp)/include" 
#LDFLAGS = -L"$(shell brew --prefix libomp)/lib"

CXX = /usr/bin/clang++
CXXFLAGS = -O3 -std=c++17 -Iinclude -Iinclude/MC -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
LDFLAGS = -L/opt/homebrew/opt/libomp/lib -lomp

#CXX = mpicxx
#CXXFLAGS = -O3 -std=c++17 -Iinclude -fopenmp

SRC_DIR = src
INC_DIR = include
OBJ_DIR = obj
BIN_DIR = bin
TARGET = $(BIN_DIR)/RNAref3D

# 再帰的に .cpp を収集
SRC_FILES := $(shell find $(SRC_DIR) -name "*.cpp")
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))
DEP_FILES := $(OBJ_FILES:.o=.d)

INC_DIRS := $(shell find $(INC_DIR) -type d)
CXXFLAGS += $(addprefix -I, $(INC_DIRS))

all: $(TARGET)

$(TARGET): $(OBJ_FILES)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@

# パスを維持しながら obj にビルド
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -MMD -MF $(OBJ_DIR)/$*.d -c $< -o $@

-include $(DEP_FILES)

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

.PHONY: all clean

