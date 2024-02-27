include ../Makefile.config

MAIN = main

.PHONY: dirs all

all: dirs main

%.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LDFLAGS)

main: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(BUILD_PATH)/$@ $^ $(LDFLAGS)

dirs:
	@echo "Creating directories"
	@mkdir -p $(BUILD_PATH)
