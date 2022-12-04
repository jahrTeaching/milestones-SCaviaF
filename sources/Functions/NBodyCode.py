from numpy import reshape, zeros
from numpy.linalg import norm

def NBody(U, t, Nb, Nc): 
    F =  zeros(len(U))
    Us  = reshape(U,(Nb, Nc, 2))
    dUs = reshape(F,(Nb, Nc, 2))

    r = reshape(Us[:,:,0],(Nb, Nc))
    v = reshape(Us[:,:,1],(Nb, Nc))
    dr = reshape(dUs[:,:,0], (Nb,Nc))
    dv = reshape(dUs[:,:,1], (Nb,Nc))
    
    for i in range(Nb):
        dr[i,:] = v[i,:]
        for j in range(Nb):
            if j != i:
                dv[i,:] = dv[i,:] + (r[j,:]-r[i,:])/norm(r[j,:]-r[i,:])**3 
    
    return F




