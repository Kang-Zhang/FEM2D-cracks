SHELL = /bin/sh
CXX = g++
INCLUDES = -I.
#CXXFLAGS = -Wall -Woverloaded-virtual -g -O0 ${INCLUDES}
CXXFLAGS = -Wall -Woverloaded-virtual -O3 ${INCLUDES}
LDFLAGS =
LD = ${CXX} ${LDFLAGS}
GLUTLIBS = -lglut -lGLU -lGL
#GLUTLIBS = -lglut -lGL -lGLU -L/usr/X11R6/lib -lXmu -lXi -lX11
#GLUTLIBS = -lglut32 -lopengl32 -lglu32
CXXSOURCES = OGLUI.cc \
	Color.cc \
	Serializable.cc \
	Font.cc \
	ToggleButtonWidget.cc \
	PushButtonWidget.cc \
	LabelWidget.cc \
	BevelWidget.cc \
	Widget.cc \
	SliderWidget.cc \
	WindowWidget.cc \
	UIWidget.cc \
	Position.cc \
	SmartLayoutManager.cc \
	Position.cc \
	Geometry.cc \
	Size.cc

OBJECTS = ${CXXSOURCES:.cc=.o}

LIBRARY = OGLUI.a

.cc.o:
	${CXX} ${CXXFLAGS} -c $<

.PHONY: all
all: depend library tests

.PHONY: library
library: ${LIBRARY}

${LIBRARY}: ${OBJECTS}
	ar cr ${LIBRARY} ${OBJECTS}

.PHONY: build
build: clean all

.PHONY:	tests
tests:	test test2 test3 test4 test5 test6 test7

test:	test.o ${LIBRARY}
	${LD} ${CXXFLAGS} -o test test.o ${LIBRARY} ${GLUTLIBS} -lm

test2:	test2.o ${LIBRARY}
	${LD} ${CXXFLAGS} -o test2 test2.o ${LIBRARY} ${GLUTLIBS} -lm

test3:	test3.o ${LIBRARY}
	${LD} ${CXXFLAGS} -o test3 test3.o ${LIBRARY} ${GLUTLIBS} -lm

test4:	test4.o ${LIBRARY}
	${LD} ${CXXFLAGS} -o test4 test4.o ${LIBRARY} ${GLUTLIBS} -lm

test5:	test5.o ${LIBRARY}
	${LD} ${CXXFLAGS} -o test5 test5.o ${LIBRARY} ${GLUTLIBS} -lm

test6:	test6.o ${LIBRARY}
	${LD} ${CXXFLAGS} -o test6 test6.o ${LIBRARY} ${GLUTLIBS} -lm

test7:	test7.o ${LIBRARY}
	${LD} ${CXXFLAGS} -o test7 test7.o ${LIBRARY} ${GLUTLIBS} -lm

.depend depend: Makefile
	@echo '...making dependencies'
	-@g++ -I. -M *.cc > .depend ${CXXFLAGS} 2> /dev/null

.PHONY: clean
clean:
	-@/bin/rm *.o
	-@/bin/rm *~
	-@/bin/rm main
	-@/bin/rm core
	-@/bin/rm core.*
	-@/bin/rm .depend
	-@/bin/rm *.exe
	-@/bin/rm ${LIBRARY}
	-@/bin/rm test
	-@/bin/rm test2
	-@/bin/rm test3
	-@/bin/rm test4
	-@/bin/rm test5
	-@/bin/rm test6 test7

sinclude .depend
