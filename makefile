TARGET=3dworld
BUILD=obj
VPATH=$(BUILD) src src/texture_tile_blend

TARGA=Targa
GLI=dependencies/gli
GLM=dependencies/glm
INCLUDES=-Isrc -Isrc/texture_tile_blend -I$(TARGA) -I$(GLI) -I$(GLM) -Idependencies/meshoptimizer/src -Idependencies/stb
DEFINES=-DENABLE_PNG -DENABLE_TIFF -DENABLE_DDS -DENABLE_STB_IMAGE -DENABLE_ASSIMP
# Note: extra warnings can be useful, but GLI and Targa generate too many warnings
CXXFLAGS=-g -Wall -O3 -fopenmp $(INCLUDES) $(DEFINES) -Wextra -Wno-unused-parameter -Wno-implicit-fallthrough \
#-Wstrict-aliasing=2 -Wunreachable-code -Wcast-align -Wcast-qual -Wsign-compare -Wsign-promo -Wdisabled-optimization -Winit-self -Wlogical-op -Wmissing-include-dirs -Wnoexcept -Woverloaded-virtual -Wredundant-decls -Wstrict-null-sentinel -Wno-unused -Wno-variadic-macros -Wno-parentheses -fdiagnostics-show-option -fasynchronous-unwind-tables -fexceptions -Werror=implicit-function-declaration -pedantic -pedantic-errors -Wformat=2 -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wimport -Winvalid-pch -Wmissing-field-initializers -Wmissing-format-attribute -Wpacked -Wpointer-arith -Wstack-protector -fstack-protector-strong -D_FORTIFY_SOURCE=2 -Wunused -Wvariadic-macros -Wwrite-strings -Werror=return-type -D_GLIBCXX_ASSERTIONS -fexceptions -fasynchronous-unwind-tables -Wctor-dtor-privacy -Wnon-virtual-dtor
OBJS=$(shell cat obj_list)

LINK=$(CPP) -fopenmp $(INCLUDES)
LDFLAGS=-lpthread `pkg-config --libs zlib libpng libtiff-4 xrender glew freealut` -lglut -fopenmp -lassimp

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
