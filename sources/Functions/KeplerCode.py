from numpy import array

"Function"

def Kepler(U, t):

    return  array([U[2],
                   U[3],
                   -U[0]/((U[0]**2 + U[1]**2)**1.5),
                   -U[1]/((U[0]**2 + U[1]**2)**1.5) ] ) 
