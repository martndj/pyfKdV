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
baseLib=spectralLib.o
kdvLib=kdvPropLib.o
wrapLib=kdvWrapC.o

kdv: $(baseLib) $(kdvLib)

#pyKdV:

%.o: %.f95
	$(compiler) -o $@ -c $< $(optflag)

clean:
	rm -rf *.o *.mod build fKdV.c

wipe:
	rm -rf *.o *.mod *.x build fKdV.c *.so
