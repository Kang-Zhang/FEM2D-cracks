SHELL = /bin/sh
CXX = g++
#OGLUIDIR = ${HOME}/Work/OGLUI
OGLUIDIR = ../OGLUI
OGLUIINC = ${OGLUIDIR}
OGLUILIBS = ${OGLUIDIR}/OGLUI.a
INCLUDES = -I. -I${OGLUIINC}
#CXXFLAGS = -Wall -g -O0 ${INCLUDES}
#CXXFLAGS = -Wall -O2 ${INCLUDES}
CXXFLAGS = -Wall -O3 ${INCLUDES}
#CXXFLAGS = -Wall -O6 ${INCLUDES}
#CXXFLAGS = -Wall -pg -O0 ${INCLUDES}
#CXXFLAGS = -Wall -g -O6 ${INCLUDES}
#CXXFLAGS = -w0 -x cxx -g ${INCLUDES}
#CXXFLAGS = -O6 -x cxx -g ${INCLUDES}
LDFLAGS =
#LDFLAGS = -pg
LD = ${CXX} ${LDFLAGS}
GLUTLIBS = -lglut -lGL -lGLU
#GLUTLIBS = -lglut32 -lopengl32 -lglu32

.cc.o:
	${CXX} ${CXXFLAGS} -c $<

.c.o:
	${CXX} ${CXXFLAGS} -c $<

OBJECTS=jo2fem.o Tokenizer.o die.o gauss_jordan.o \
	Matrix.o pinv.o TriangulatorDelaunay2D.o Vector3d.o \
	StressTensor.o realcubicroots.o nullspace.o xmath.o

.PHONY: all
all: depend jo2fem

.PHONY: build
build: clean all

jo2fem:	${OBJECTS} ${OGLUILIBS}
	${CXX} -o jo2fem ${OBJECTS} ${OGLUILIBS} ${GLUTLIBS} -lm

.PHONY: depend
.depend depend: Makefile
	@echo '...making dependencies'
	-@g++ -I. -M *.cc *.c > .depend ${CXXFLAGS} 2> /dev/null

clean:
	-@/bin/rm *.o
	-@/bin/rm *~
	-@/bin/rm core
	-@/bin/rm jo2fem
	-@/bin/rm .depend
	-@/bin/rm *.exe

sinclude .depend
