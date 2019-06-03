CPP=g++

BUILD=obj
TARGA=Targa
GLI=dependencies/gli
CPPFLAGS=-g -Wall -O3 -fopenmp -Isrc -Isrc/texture_tile_blend -I$(TARGA) -I$(GLI) -DENABLE_JPEG -DENABLE_PNG -DENABLE_TIFF -DENABLE_DDS
TARGET=3dworld
OBJS=$(shell cat obj_list)
VPATH=$(BUILD):src:src/texture_tile_blend

LINK=g++ $(CPPFLAGS)
LFLAGS=-lz -lpng -ljpeg -ltiff -lpthread $(shell pkg-config --libs xrender) -lglut -lGLEW -lGLU -lGL -lopenal -lalut

#  In most cases, you should not change anything below this line.
all : $(TARGET)

# disable old-style .SUFFIXES rules. this may not be needed?
.SUFFIXES:

%.o : %.C $(BUILD)/%.d
	$(CPP) $(CPPFLAGS) -MMD -c $< -o $(abspath $(BUILD)/$@)
%.o : %.cc $(BUILD)/%.d
	$(CPP) $(CPPFLAGS) -MMD -c $< -o $(abspath $(BUILD)/$@)
%.o : %.cpp $(BUILD)/%.d
	$(CPP) $(CPPFLAGS) -MMD -c $< -o $(abspath $(BUILD)/$@)

# DEPENDENCIES = $(OBJS:.o=.d)

# 
# Targets:
# 

$(TARGET): $(OBJS)
	cd $(BUILD) && $(LINK) $(FLAGS) -o $(TARGET) $(OBJS) $(LFLAGS)

.PHONY : clean
clean:
	-rm -fr $(BUILD)

$(OBJS): | $(BUILD)
$(BUILD):
	@mkdir -p $(BUILD)

$(BUILD)/%.d: ;
.PRECIOUS: $(BUILD)/%.d
-include $(patsubst %,$(BUILD)/%.d,$(basename $(OBJS)))
