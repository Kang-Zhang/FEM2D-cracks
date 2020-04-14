SHELL = /bin/sh

all:
	+@(cd OGLUI ; make )
	+@(cd Leaf ; make )
	+@(cd Relax ; make )
	+@(cd 2D ; make )

clean:
	+@(cd OGLUI ; make clean )
	+@(cd Leaf ; make clean )
	+@(cd Relax ; make clean )
	+@(cd 2D ; make clean )
	+/bin/rm *~
