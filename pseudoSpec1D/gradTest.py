import numpy as np

def gradientTest(x, fct, grad, powRange=[-1,-14],
                output=True, args=()):
    if (not callable(fct) or not callable(grad)):
        raise TypeError()
    J0=fct(x, *args)
    gradJ0=grad(x,  *args)
    n2GradJ0=np.dot(gradJ0, gradJ0)

    test={}
    for power in xrange(powRange[0],powRange[1], -1):
        eps=10.**(power)
        Jeps=fct(x-eps*gradJ0, *args)
        res=((J0-Jeps)/(eps*n2GradJ0))
        test[power]=[Jeps, res]

    if output:  print(gradTestString(J0, n2GradJ0, test))
    return (J0, n2GradJ0, test)
    
def gradTestString(J0, n2GradJ0, test):
    s="----| Gradient test |------------------\n"
    s+="  J0      =%+25.15e\n"%J0
    s+=" |grad|^2 =%+25.15e\n"%n2GradJ0
    for i in  (np.sort(test.keys())[::-1]):
        s+="%4d %+25.15e  %+25.15e\n"%(i, test[i][0], test[i][1])
    return s

#--------------------------------------------------------------------
#====================================================================
#--------------------------------------------------------------------
