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

optflag=-fbounds-check -fPIC -O3
libflag=-L$(libDir) -larpack -lfftw3 -lm
baseLib=matrixLib.o spectralLib.o
kdvLib=kdvPropLib.o kdvPropLib_pt.o kdvTLMPropLib.o kdvTLMPropLib_pt.o kdvLanczosLib.o kdvTestLib.o 
kdvTests=kdvTestAdj.x kdvTestGrad.x  
wrapLib=kdvWrapC.o 


kdv: $(baseLib) $(kdvLib) 

testGrad: $(baseLib) $(kdvLib) $(kdvTests)
	./kdvTestGrad.x

testAdj: $(baseLib) $(kdvLib) $(kdvTests)
	./kdvTestAdj.x

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
