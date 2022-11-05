from numpy import array

"Linear Oscilator function"

def LinOsc(U, t):

    return array([U[1], -U[0]])
