TARGET=3dworld
BUILD=obj
VPATH=$(BUILD) src src/texture_tile_blend

TARGA=Targa
GLI=dependencies/gli
INCLUDES=-Isrc -Isrc/texture_tile_blend -I$(TARGA) -I$(GLI) -I$(GLI)/external
DEFINES=-DENABLE_JPEG -DENABLE_PNG -DENABLE_TIFF -DENABLE_DDS
# Note: extra warnings can be useful, but GLI and Targa generate too many warnings
CXXFLAGS=-g -Wall -O3 -fopenmp $(INCLUDES) $(DEFINES) #-Wextra -Wno-unused-parameter -Wno-implicit-fallthrough
OBJS=$(shell cat obj_list)

LINK=$(CPP) -fopenmp $(INCLUDES)
LDFLAGS=-lpthread `pkg-config --libs zlib libpng libjpeg libtiff-4 xrender glew freealut` -lglut -fopenmp

# For creating dependencies files
DEPFLAGS = -MT $@ -MMD -MP -MF $(BUILD)/$*.Td
POSTCOMPILE = mv -f $(BUILD)/$*.Td $(BUILD)/$*.d

# Change the verbosity of the makefile.
ifeq ($(VERBOSE),1)
Q :=
else
Q := @
endif

# Compile 3dworld
all: $(TARGET)

# Link the target
$(TARGET): $(OBJS)
	@echo "Linking $<"
	$(Q)cd $(BUILD) && $(CXX) $(INCLUDES) -o $(TARGET) $(OBJS) $(LDFLAGS)

# Compile source files
%.o : %.cpp $(BUILD)/%.d
	@echo "Compiling $<"
	$(Q)$(CXX) $(DEPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c $(abspath $<) -o $(abspath $(BUILD)/$@)
	@$(POSTCOMPILE)

# Delete compiled files
.PHONY: clean
clean:
	-rm -fr $(BUILD)

# Create the directory before compiling sources
$(OBJS): | $(BUILD)
$(BUILD):
	@mkdir -p $(BUILD)

# Create the dependency files
$(BUILD)/%.d: ;
.PRECIOUS: $(BUILD)/%.d

-include $(patsubst %,$(BUILD)/%.d,$(basename $(OBJS)))
