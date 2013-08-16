...In translation from fKdV4
============================

pyfKdV is a fortran 95 propagator wrapped in python through cython compiler.

__In active development__:

 * TLM nor adjoint implemented in this version
 * nor Lanczos, singular vector calculation
 * fKdV4 is my stable version for now. It is an (almost) Object-orientated Fortran 95 integrator for the Korteweg de-Vries differential equation.  It also contain a (kind of) wrapper in python, but it is very clumsy (passing everything by file transfer between python and fortran...)

It aims at calculating singular vectors using Lanczos algorithm.

About augmented KdV system and it's use in atmospheric dynamics:

 * [Patoine, A., T. Warn, 1982: The Interaction of Long, Quasi-Stationary Baroclinic Waves with Topography. J. Atmos. Sci., 39, 1018–1025.](http://journals.ametsoc.org/doi/abs/10.1175/1520-0469(1982\)039%3C1018%3ATIOLQS%3E2.0.CO%3B2)
 * [Hodyss, D. and Nathan, T. R. Solitary rossby waves in zonally varying jet ows. Geophys. Astrophys., 262, 2002.](http://www.tandfonline.com/doi/abs/10.1080/03091920290011012#.Ug1egSHPTMU)

