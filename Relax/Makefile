SHELL = /bin/sh
CXX = g++
INCLUDES = -I.
CXXFLAGS = -Wall -O3 ${INCLUDES}
LDFLAGS =
LD = ${CXX} ${LDFLAGS}
GLUTLIBS = -L/usr/X11R6/lib -lglut -lGL -lGLU
#GLUTLIBS = -lglut32 -lglu32 -lopengl32

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

.cc.o:	Makefile
	${CXX} ${CXXFLAGS} -c $<

.PHONY: all
all: relax relax-anim

.PHONY: nothing

.PHONY: build
build:
	make clean
	make all

relax: relax.o
	${LD} -o relax $^ ${GLUTLIBS} -lm

relax-anim: relax-anim.o
	${LD} -o relax-anim $^ ${GLUTLIBS} -lm -lpng

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
	-@/bin/rm -f core
	-@/bin/rm -f relax
	-@/bin/rm -f relax-anim
	-@/bin/rm -f *.exe
	-@/bin/rm -f *.d
	-@/bin/rm -f *.stackdump
	-@/bin/rm allsources.txt
