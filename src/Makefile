#! /bin/sh

#--------------------------------------------------------------------#
#  Makefile pour cKdV5
#
#	makefile references
# 	http://gl.developpez.com/tutoriel/outil/makefile/
#
#--------------------------------------------------------------------#


#----| variables |---------------------------------------------------#
user=martn
libDir=/usr/lib
pyKdVDir=../pyKdV/
userLibDir=/home/${user}/lib/
compiler=gfortran

optflag=-fbounds-check -fPIC
libflag=-L$(libDir) -larpack -lfftw3 -lm
baseLib=matrixLib.o spectralLib.o
kdvLib=kdvPropLib.o kdvTLMPropLib.o kdvTLMTestLib.o kdvLanczosLib.o 
kdvTests=kdvTestAdj.x kdvTestGrad.x 
wrapLib=kdvWrapC.o 

kdv: $(baseLib) $(kdvLib)

tests:  $(baseLib) $(kdvLib) $(kdvTests)
	./kdvTestAdj.x
	./kdvTestGrad.x

pyKdV: libKdV.so
	python setup.py build_ext --inplace
	cp fKdV.so $(pyKdVDir)
	cp libKdV.so $(userLibDir)
    

%.o: %.f95
	$(compiler) -o $@ -c $< $(optflag)

%.so: $(baseLib) $(kdvLib) $(wrapLib)
	$(compiler) -shared -o $@ $? $(libflag) $(optflag)

%.x: $(baseLib) $(kdvLib)
	$(compiler) -c $*.f95 $(optflag)
	$(compiler) $*.o $? -o $@ $(libflag) $(optflag)
	chmod +x $@

clean:
	rm -rf *.o *.mod build fKdV.c

wipe:
	rm -rf *.o *.mod *.x build fKdV.c *.so *.pyc
