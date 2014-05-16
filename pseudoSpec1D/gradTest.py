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

def gradTestCheck(gradTest, verbose=False):
    """
    Inspect a gradient test result
    and return the number of significant digit
        
        consider cases of 0.999... and of 1.000...
    """
    def significantFigures(testValues):
        def countFigures(figure, digits):
            digits=digits.split('.')[1]
            pre=True
            count=0
            n=0
            while n < len(digits) and pre==True:
                if int(digits[n])==figure:
                    count+=1
                else:
                    pre=False
                n+=1
            return count

        significant={}
        for p in sorted(testValues.keys())[::-1]:
            significant[p]=0
            digits=str(testValues[p][1])
            if digits=='nan':
                pass
            elif int(digits.split('.')[0])==0:
                significant[p]=countFigures(9, digits)
            elif int(digits.split('.')[0])==1:
                significant[p]=countFigures(0, digits)
            else:
                pass
        return significant
    
    significant=significantFigures(gradTest[2])
    valid=0
    count=0
    i=0
    for p in  sorted(significant.keys())[::-1]:
        if significant[p]>count: 
            count+=1
        elif significant[p]==count:
            pass
        elif significant[p]==0:
            count=0
        else:
            count=1
        if count>valid: valid=count
        if verbose :  print(significant[p], count, valid)
    return valid, significant
#--------------------------------------------------------------------
#====================================================================
#--------------------------------------------------------------------
