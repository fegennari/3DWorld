VPATH=../src
GLUI=../../glui
TARGA=../Targa
GLI=../../gli
CPPFLAGS=-g -Wall -O3 -fopenmp -I$(TARGA) -I$(GLUI)/include -I$(GLI) -I../src/texture_tile_blend -DENABLE_JPEG -DENABLE_PNG -DENABLE_TIFF -DENABLE_DDS
TARGET=../obj/3dworld
OBJS=$(shell cat ../obj_list)

CPP=g++
#LINK=$(CPP) $(CPPFLAGS) -lz -lpng -lpthread -L/usr/X11R6/lib64 -lglut -lGL -lGLU
LINK=$(CPP) $(CPPFLAGS) -L$(GLUI)/lib

LFLAGS=-lz -lpng -ljpeg -ltiff -lpthread $(shell pkg-config --libs xrender) -lglui -lglut -lGLEW -lGLU -lGL -lopenal -lalut

#  In most cases, you should not change anything below this line.
ifeq ($(shell test -L makefile ; echo $$? ),1)
all : 
	@echo "makefile should be a symbolic link to avoid accidentally building in the src directory ... attempting to create ../obj,../lib,../run, symlink makefile in ../obj, and recurse make into ../obj"
	-mkdir ../obj
	-mkdir ../run
	-ln -s ../makefile ../obj/makefile
	cd ../obj && $(MAKE)
else

# disable old-style .SUFFIXES rules. this may not be needed?
.SUFFIXES:

%.o : %.C
	$(CPP) $(CPPFLAGS) -MMD -c $<
%.o : %.cc
	$(CPP) $(CPPFLAGS) -MMD -c $<
%.o : %.cpp
	$(CPP) $(CPPFLAGS) -MMD -c $<

%.d: %.C
	touch $@
%.d: %.cc
	touch $@
%.d: %.cpp
	touch $@

DEPENDENCIES = $(OBJS:.o=.d)

# 
# Targets:
# 

all : $(TARGET)

$(TARGET): $(OBJS)
	$(LINK) $(FLAGS) -o $(TARGET) $(OBJS) $(LFLAGS)

.PHONY : clean

clean:
	-rm -f $(TARGET) $(OBJS) $(DEPENDENCIES) make.dep

make.dep: $(DEPENDENCIES)
	-cat $(DEPENDENCIES) > make.dep

include make.dep

endif
