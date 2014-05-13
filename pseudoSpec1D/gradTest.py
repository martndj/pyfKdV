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
if __name__=='__main__':

    import pyKdV as kdv
    
    grid=kdv.PeriodicGrid(144)
    param=kdv.Param(grid)
    dt=0.01
    modelNL=kdv.kdvLauncher(param, dt=dt)
    tlm=kdv.kdvTLMLauncher(param)
    
    nDt=500
    tInt=nDt*dt
    u0=kdv.rndSpecVec(grid, Ntrc=12, seed=1)
    u=modelNL.integrate(u0, tInt)
    mTLM=kdv.kdvTLMLauncher(modelNL.param)
    mTLM.reference(u)
    freq=3
    nDtList=[i*nDt/freq for i in xrange(1,freq+1)]
    tList=dt*np.array(nDtList)

    def fct1(u0, times, model, tlm):
        d_u=model.d_intTimes(u0, times)
        J=0.
        for t in d_u.keys():
            J+=0.5*np.dot(d_u[t], d_u[t])
        return J
    def gradF1(u0, times, model, tlm):
        u=model.integrate(u0, times[-1])
        tlm.reference(u)
        d_u=model.d_intTimes(u0, times)
        return tlm.d_intTimesAdj(d_u)
    
    kdv.gradientTest(u0, fct1, gradF1, args=(tList, modelNL, mTLM))
    
    def fct2(u0, nDtL, model, tlm):
        d_u=model.d_nDtInt(u0, nDtL)
        J=0.
        for t in d_u.keys():
            J+=0.5*np.dot(d_u[t], d_u[t])
        return J
    def gradF2(u0, nDtL, model, tlm):
        #u=model.integrate(u0, nDtL[-1]*model.dt)
        #tlm=kdv.kdvTLMLauncher(model.param, traj=u)
        d_u=model.d_nDtInt(u0, nDtL)
        return tlm.d_nDtIntAdj(d_u)
    
    kdv.gradientTest(u0, fct2, gradF2, args=(nDtList, modelNL, mTLM))
    
