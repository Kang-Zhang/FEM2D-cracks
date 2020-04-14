SHELL = /bin/sh
CXX = g++
GSLINCLUDES = `gsl-config --cflags`
GSLLIBS = `gsl-config --libs`
#CXX = cxx
#INCLUDES = ${GSLINCLUDES}
INCLUDES = -I.
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
#GLUTLIBS = -L/home/users/federl/usr/lib -L/usr/X11R6/lib -lglut -lGL -lGLU
#GLUTLIBS = -L/usr/X11R6/lib -lglut -lGL -lGLU -lXi -lXmu
#GLUTLIBS = -L/usr/lib -lglut -lGL -lGLU -L /usr/X11R6/lib -lXmu -lXi
#GLUTLIBS = -lglut32 -lglu32 -lopengl32
GLUTLIBS = -lglut -lGLU -lGL
MAINSOURCES = main.cc Matrix.cc \
	  conjugate_gradient.cc \
	  SparseMatrix.cc \
	  Tokenizer.cc \
	  Distribution.cc \
	  Map.cc \
	  Image.cc \
	  plane_split.cc \
	  StressTensor.cc \
	  nullspace.cc \
	  xmath.cc \
	  Vector3d.cc \
	  projected_area.cc \
	  realcubicroots.cc \
	  gauss_jordan.cc \
	  pinv.cc \
	  die.cc \
	  readrgb.cc


.PHONY: default
default: all

.PHONY:
allsources.txt:
	echo 'ALLSOURCES = ' *.cc > allsources.txt

sinclude allsources.txt

sinclude ${ALLSOURCES:.cc=.d}

%.d: %.cc
	-@echo -n 'Making dependency file for ' $< '...'
	-@set -e; $(CXX) -MM -nostdinc $(CXXFLAGS) $< \
		| sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; \
		[ -s $@ ] || rm -f $@
	-@echo ' ok'

OTHERSOURCES = viewer.cc viewer-cyl.cc viewsurface.cc

MAINOBJECTS = ${MAINSOURCES:.cc=.o}

OTHEROBJECTS = ${OTHERSOURCES:.cc=.o}

.cc.o:	Makefile
	${CXX} ${CXXFLAGS} -c $<

.PHONY: all
all: main viewer viewer-cyl genplane genplane2 genplane3\
	relax relax-anim gencylinder gensphere perlin \
	surface2fem gencylinder2 genplane-bumps

.PHONY: nothing

.PHONY: build
build:
	make clean
	make all

main:	${MAINOBJECTS}
	${LD} ${CXXFLAGS} -o main $^ ${GLUTLIBS} -lpng -lm

jo2fem:	jo2fem.o jo2fem.o Tokenizer.o die.o gauss_jordan.o \
	Matrix.o pinv.o TriangulatorDelaunay2D.o
	${CXX} -o jo2fem $^ ${GLUTLIBS} -lm

viewer:	viewer.o die.o Tokenizer.o Vector3d.o gl2ps.o
	${LD} -o viewer $^ ${GLUTLIBS} -lm

viewer-cyl: viewer-cyl.o die.o Tokenizer.o
	${LD} -o viewer-cyl $^ ${GLUTLIBS} -lm

viewsurface: viewsurface.o die.o
	${LD} -o viewsurface $^ ${GLUTLIBS} -lm

genplane: die.o genplane.o Distribution.o 
	${LD} -o genplane $^ -lm

genplane2: die.o genplane2.o Distribution.o MapTexture2D.o Image.o readrgb.o
	${LD} -o genplane2 $^ -lm

genplane3: die.o genplane3.o Distribution.o MapTexture2D.o Image.o readrgb.o
	${LD} -o genplane3 $^ -lm

genplane-bumps: die.o genplane-bumps.o Distribution.o 
	${LD} -o genplane-bumps $^ -lm

gencylinder: gencylinder.o
	${LD} -o gencylinder $^ -lm

gencylinder2: gencylinder2.o Distribution.o die.o
	${LD} -o gencylinder2 $^ -lm

gensphere: gensphere.o
	${LD} -o gensphere $^ -lm

relax: relax.o
	${LD} -o relax $^ ${GLUTLIBS} -lm

relax-anim: relax-anim.o
	${LD} -o relax-anim $^ ${GLUTLIBS} -lm -lpng

perlin:	perlin.o
	${LD} -o perlin $^ -lm

eigentest: eigentest.o eigen.o Matrix.o die.o
	${CXX} -o eigentest $^ -lm

stressTensorTest: stressTensorTest.o StressTensor.o nullspace.o xmath.o
	${CXX} -o stressTensorTest $^ -lm

planetest: plane_split.o
	${CXX} -o planetest $^ -DSELFTEST -lm

surface2fem:	surface2fem.o Tokenizer.o die.o
	${CXX} -o surface2fem $^ -lm

.PHONY: removedeps
removedeps:
	-@/bin/rm -f *.d

.PHONY: makedeps
makedeps: ${ALLSOURCES:.cc=.d}

.PHONY: depend
depend:
	make removedeps
	make makedeps

.PHONY: clean
clean:
	-@/bin/rm -f *.o
	-@/bin/rm -f *~
	-@/bin/rm -f main
	-@/bin/rm -f core
	-@/bin/rm -f genplane
	-@/bin/rm -f genplane2
	-@/bin/rm -f genplane3
	-@/bin/rm -f genplane-bumps
	-@/bin/rm -f gencylinder
	-@/bin/rm -f gencylinder2
	-@/bin/rm -f gensphere
	-@/bin/rm -f relax
	-@/bin/rm -f relax-anim
	-@/bin/rm -f perlin
	-@/bin/rm -f viewer
	-@/bin/rm -f viewer-cyl
	-@/bin/rm -f viewsurface
	-@/bin/rm -f eigentest
	-@/bin/rm -f surface2fem
	-@/bin/rm -f jo2fem
	-@/bin/rm -f *.exe
	-@/bin/rm -f *.d
	-@/bin/rm -f *.stackdump
	-@/bin/rm allsources.txt
