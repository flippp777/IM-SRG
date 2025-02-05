# Compiler und Flags
CXX = g++
CXXFLAGS = -fdiagnostics-color=always
LDFLAGS =  -fopenmp

# Verzeichnisse
BUILD_DIR = build
BIN_DIR = $(BUILD_DIR)/bin

# Alle .cpp-Dateien sammeln
SRCS = $(wildcard src/*.cpp)
OBJS = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(SRCS))
TARGET = $(BIN_DIR)/im_srg

# Automatische Abhängigkeitsdateien generieren
DEPFLAGS = -MMD -MP
CXXFLAGS += $(DEPFLAGS)

# Alle Abhängigkeitsdateien einbinden
-include $(OBJS:.o=.d)

# Standard-Ziel: Erstelle das Programm
all: $(TARGET)

# Linken der .o-Dateien zum fertigen Programm
$(TARGET): $(OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET) $(LDFLAGS)

# Kompilieren jeder .cpp-Datei zu einer .o-Datei
$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Stelle sicher, dass die Verzeichnisse existieren
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Clean-Target, um die Build-Dateien zu entfernen
clean:
	rm -rf $(BUILD_DIR)

