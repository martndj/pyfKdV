import numpy as np

from pseudoSpec1D import PeriodicGrid, Trajectory,\
                            Launcher, TLMLauncher

class ValidateTLMError(Exception):
    pass


def gradTest(tlm, nlModel, ic, tInt,
                output=True, powRange=[-1,-14]):

    if not (isinstance(tlm,TLMLauncher)):
        raise ValidateTLMError("tlm <TLMLauncher>")
    if not (isinstance(nlModel,Launcher)):
        raise ValidateTLMError("nlModel <Launcher>")
    if not tlm.isInitialized:
        raise ValidateTLMError(
                    "Not initialized with a reference trajectory")

    def evolvedNorm2(x):
        xf=nlModel.integrate(x, tInt).final
        return 0.5*np.dot(xf,xf)
    #  grad:
    #       tlm.adjoint(nl.integrate(x, tInt), tInt)


    traj=nlModel.integrate(ic, tInt)
    final=traj.final
    J0=evolvedNorm2(final)

    tlm.initialize(traj)
    gradJ0=tlm.adjoint(final, tInt)
    n2GradJ0=np.dot(gradJ0, gradJ0)

    test={}
    for power in xrange(powRange[0],powRange[1], -1):
        eps=10.**(power)
        Jeps=evolvedNorm2(nlModel.integrate(ic-eps*gradJ0, tInt).final)
        
        res=((J0-Jeps)/(eps*n2GradJ0))
        test[power]=[Jeps, res]

    if output:
        print("----| Gradient test |------------------")
        print("  M(ic)      =%+25.15f"%J0)
        print(" |grad|^2 =%+25.15f"%n2GradJ0)
        for i in  (np.sort(test.keys())[::-1]):
            print("%4d %+25.15f  %+25.15f"%(i, test[i][0], test[i][1]))

    return (J0, n2GradJ0, test)

