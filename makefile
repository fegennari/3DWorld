
VPATH=../src
#CPPFLAGS=-g -Wall 
CPPFLAGS=-g -Wall -O3
TARGET=../lib/3dworld
OBJS=$(shell cat ../src/obj_list)
TARGET2=
OBJS2= 

LIB_TARGET=
LIB_OBJS=

#LINK=g++ $(CPPFLAGS) -lz -lpng -lpthread -L/usr/X11R6/lib64 -lglut -lGL -lGLU
LINK=g++ $(CPPFLAGS) -lz -lpng -lpthread $(shell pkg-config --libs xrender) -lglui -lglut -lGL -lGLU

LFLAGS=$(LIB_TARGET)

#  In most cases, you should not change anything below this line.
ifeq ($(shell test -L makefile ; echo $$? ),1)
all : 
	@echo "makefile should be a symbolic link to avoid accidentally building in the src directory ... attempting to create ../obj,../lib,../run, symlink makefile in ../obj, and recurse make into ../obj"
	-mkdir ../obj
	-mkdir ../run
	-mkdir ../lib
	-ln -s ../src/makefile ../obj/makefile
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

CPP=g++

%.d: %.C
	touch $@
%.d: %.cc
	touch $@
%.d: %.cpp
	touch $@

DEPENDENCIES = $(OBJS2:.o=.d) $(OBJS:.o=.d) $(LIB_OBJS:.o=.d)

# 
# Targets:
# 

all : $(TARGET) $(TARGET2)

ifneq ($(LIB_TARGET),"")
$(LIB_TARGET): $(LIB_OBJS) 
	$(LINK) $(FLAGS) -shared -o $(LIB_TARGET) $(LIB_OBJS)
endif

$(TARGET): $(OBJS) $(LIB_TARGET)
	$(LINK) $(FLAGS) -o $(TARGET) $(OBJS) $(LFLAGS)

ifneq ($(TARGET2),"")
$(TARGET2): $(OBJS2) $(LIB_TARGET)
	$(LINK) $(FLAGS) -o $(TARGET2) $(OBJS2) $(LFLAGS)
endif

.PHONY : clean

clean:
	-rm -f $(TARGET) $(TARGET2) $(OBJS) $(OBJS2) $(LIB_TARGET) $(LIB_OBJS) $(DEPENDENCIES) make.dep

make.dep: $(DEPENDENCIES)
	-cat $(DEPENDENCIES) > make.dep

include make.dep

endif