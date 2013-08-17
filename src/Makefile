#! /bin/sh

#--------------------------------------------------------------------#
#  Makefile pour cKdV5
#
#	makefile references
# 	http://gl.developpez.com/tutoriel/outil/makefile/
#
#--------------------------------------------------------------------#


#----| variables |---------------------------------------------------#
libDir=/usr/lib
compiler=gfortran

optflag=-fbounds-check -fPIC
libflag=-L$(libDir) -larpack -lfftw3 -lm
baseLib=matrixLib.o spectralLib.o
kdvLib=kdvPropLib.o kdvTLMPropLib.o kdvTLMTestLib.o 
kdvTests=kdvTestAdj.x
wrapLib=kdvWrapC.o

kdv: $(baseLib) $(kdvLib)

tests:  $(baseLib) $(kdvLib) $(kdvTests)
	./kdvTestAdj.x

pyKdV: kdvLib.so
	python setup.py build_ext --inplace

%.o: %.f95
	$(compiler) -o $@ -c $< $(optflag)

%.so: $(baseLib) $(kdvLib) $(wrapLib)
	$(compiler) -shared -o $@ $? $(optflag)

%.x: $(baseLib) $(kdvLib)
	$(compiler) -c $*.f95 $(optflag)
	$(compiler) $*.o $? -o $@ $(libflag) $(optflag)
	chmod +x $@

clean:
	rm -rf *.o *.mod build fKdV.c

wipe:
	rm -rf *.o *.mod *.x build fKdV.c *.so *.pyc
