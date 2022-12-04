from ODEs.SchemesCode import Euler, RK4, CN, EulerInv, LF
from numpy import meshgrid, array, zeros, size, absolute, sqrt

"Absolute Stability Region"

def Regions(x,y,Scheme):  
    Z = zeros([size(x),size(x)],dtype=complex)
    for i in range(size(x)):
        for j in range(size(x)):
            Z[size(x)-1-j,i] = complex(x[i],y[j])

    return absolute(array(StabPoly(Scheme, Z)))



def StabPoly(Scheme, w):
    if Scheme == Euler:
        r = 1 + w
    elif Scheme == EulerInv:
        r = 1/(1-w)
    elif Scheme == CN:
        r = (1+w/2)/(1-w/2)
    elif Scheme == RK4:
        r = 1 + w + (w**2)/2 + (w**3)/6 + (w**4)/(4*3*2)
    elif Scheme == LF:
        r = sqrt(1)

    return r



