# simulation requires a function
# which enforces boundary conditions 
# on each iteration
from numba import njit,prange
from numpy import argwhere

def flow(domain,gx,gy): # flow with obstructions
    obs_idx = argwhere(domain[0]==False)
    @njit(parallel=True)
    def _bc(Ux,Uy,p):
        for ij in prange(len(obs_idx)):
            i,j=obs_idx[ij]
            p[i+1,j+1]=p[i-1,j-1]=p[i-1,j+1]=p[i+1,j-1]
            Ux[i,j]=Uy[i,j]=0
        Ux[:,0]=Ux[:,1]
        Uy[:,0]=Uy[:,1]
        p[:,0]=p[:,1]
        Ux[:,-1]=Ux[:,-2]
        Uy[:,-1]=Uy[:,-2]
        p[:,-1]=p[:,-2]
        Ux[0,:]=Ux[1,:]
        Uy[0,:]=Uy[1,:]
        p[0,:]=p[1,:]
        Ux[-1,:]=Ux[-2,:]
        Uy[-1,:]=Uy[-2,:]
        p[-1,:]=p[-2,:]
        return Ux,Uy,p
    return _bc,gx,gy

def lid(Ux): # standard lid-driven cavity problem
    _=Ux
    def bc(Ux,Uy,p):
        p[-1,:]=p[-2,:]
        p[0,:]=p[1,:]
        p[:,0]=p[:,1]
        p[:,-1]=0
        Ux[-1,:]=Uy[-1,:]=0
        Ux[0,:]=Uy[0,:]=0
        Ux[:,0]=Uy[:,0]=0
        Ux[:,-1]=_
        Uy[:,-1]=0
        return Ux,Uy,p
    return bc,0,0

def channel(gx): # channel flow
    def bc(Ux,Uy,p):
        Ux[:,0]=0
        Uy[:,0]=0
        p[:,0]=p[:,1]=p[:,2]
        Ux[:,-1]=0
        Uy[:,-1]=0
        p[:,-1]=p[:,-2]=p[:,-3]
        Ux[0,:]=Ux[1,:]
        Uy[0,:]=Uy[1,:]
        p[0,:]=p[1,:]
        Ux[-1,:]=Ux[-2,:]
        Uy[-1,:]=Uy[-2,:]
        p[-1,:]=p[-2,:]
        return Ux,Uy,p
    return bc,gx,0