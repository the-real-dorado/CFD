from numba import njit,prange
from numpy import ones,zeros,linspace,meshgrid,stack, inf
from numpy import sum as arraysum
from numpy.linalg import norm
from matplotlib.pyplot import figure,contourf,axis,show
from time import time
from .utilities import importobj,in_cli,status

@njit(parallel=True) # numba go brrr
def _NS(Ux,Uy,p,gx,gy,nu,Nx,Ny,dx,dy,dt,bc,p_it): # solving the incompressible Navier-Stokes equations
    Ux_,Uy_=Ux.copy(),Uy.copy()                    # with a poisson equation for the pressure field
    for _ in range(p_it):                          # using the intermediate velocity
        p_=p.copy()                                # by Jacobi iteration
        for i in prange(1,Nx-1):     
            for j in prange(1,Ny-1): 
                p[i,j]=( ( (p_[i+1,j]+p_[i-1,j])*dy**2 + (p_[i,j+1]+p_[i,j-1])*dx**2 ) / ( 2*(dx**2+dy**2) )
                        - ( ( dx**2*dy**2 ) / ( 2*(dx**2+dy**2) ) )
                        * ( (1/dt)*((Ux_[i+1,j]-Ux_[i-1,j])/(2*dx)+(Uy_[i,j+1]-Uy_[i,j-1])/(2*dy))
                        - ((Ux_[i+1,j]-Ux_[i-1,j])/(2*dx))**2
                        - ((Uy_[i,j+1]-Uy_[i,j-1])/(2*dy))**2
                        -2*((Ux_[i,j+1]-Ux_[i,j-1])/(2*dy))*((Uy_[i+1,j]-Uy_[i-1,j])/(2*dx)) ) )
    for i in prange(1,Nx-1):     
        for j in prange(1,Ny-1):
            Ux[i,j]=( Ux_[i,j]
                   - Ux_[i,j] * (dt/dx) * (Ux_[i,j]-Ux_[i-1,j]) - Uy_[i,j] * (dt/dy) * (Ux_[i,j]-Ux_[i,j-1])
                   - (dt/(2*dx)) * (p[i+1,j]-p[i-1,j])
                   + nu * (dt/dx**2) * (Ux_[i+1,j]-2*Ux_[i,j]+Ux_[i-1,j]) + nu * (dt/dy**2) * (Ux_[i,j+1]-2*Ux_[i,j]+Ux_[i,j-1])
                   + dt*gx )
            Uy[i,j]=( Uy_[i,j]
                   - Ux_[i,j] * (dt/dx) * (Uy_[i,j]-Uy_[i-1,j]) - Uy_[i,j] * (dt/dy) * (Uy_[i,j]-Uy_[i,j-1])
                   - (dt/(2*dy)) * (p[i,j+1]-p[i,j-1])
                   + nu * (dt/dx**2) * (Uy_[i+1,j]-2*Uy_[i,j]+Uy_[i-1,j]) + nu * (dt/dy**2) * (Uy_[i,j+1]-2*Uy_[i,j]+Uy_[i,j-1])
                   + dt*gy )
    Ux,Uy,p = bc(Ux,Uy,p)
    return Ux,Uy,p,Ux_,Uy_

def domain(file,dl,Nx,Ny,*,display=False): # domain for fluid
    if in_cli(): print('------------------------------------')
    print('*Meshing*')
    M=zeros((Nx+1,Ny+1))
    if file: # complex domain from .obj file
        D=importobj(file,dl)
        for c in D:
            M[c[0],c[1]]=1
        M[:,0]=M[:,-1]=M[0,:]=M[-1,:]=1
    if display: # render domain
        if in_cli(): print('(close window to continue)',end='\r')
        x,y=linspace(0,Nx*dl,Nx+1),linspace(0,Ny*dl,Ny+1)
        figure(figsize=(Nx/Ny*5,5),tight_layout=True)
        contourf(x,y,-M.T,cmap='bone_r')
        axis('scaled')
        axis('off')
        show()
    if in_cli(): print('------------------------------------')
    return M,dl,Nx+1,Ny+1 # passing number of points instead of divisions

class fluid: # fluid for use in simulation
    def __init__(self,domain,*,nu=1e-6):
        self.domain,self.dl,self.Nx,self.Ny=domain
        self.nu=nu
        self.x=linspace(0,(self.Nx-1)*self.dl,self.Nx)
        self.y=linspace(0,(self.Ny-1)*self.dl,self.Ny)
        self.X,self.Y=meshgrid(self.x,self.y)
    def update(self,Ux,Uy,p):
        self.Ux=Ux
        self.Uy=Uy
        self.p=p
    def stats(self,it,re,t):
        self.it=it
        self.re=re
        self.t=t

def simulation(fluid,model,*,dt=1e-6,p_it=200,rt=5e-4,max_it=10000): # calls calculation function until convergence or max iterations
    print('*Initializing*')
    F=fluid
    Nx,Ny=F.Nx,F.Ny
    dx=dy=F.dl
    Ux,Uy=zeros((Nx,Ny)),zeros((Nx,Ny))
    p=zeros((Nx,Ny)) # why one???
    nu=F.nu
    bc,gx,gy = model
    Ux,Uy,p=bc(Ux,Uy,p)
    Ux,Uy,p,Ux_,Uy_=_NS(Ux,Uy,p,gx,gy,nu,Nx,Ny,dx,dy,dt,bc,p_it)
    re=convergence(Ux,Uy,Ux_,Uy_)
    print(f'time step:       {dt}s')
    print(f'rel. tolerance:  {rt}')
    print(f'Jacobi loops:    {p_it}')
    print(f'max. iterations: {max_it}')
    cli=in_cli()
    if cli: print('------------------------------------')
    print(f'*Solving*')
    ts=time()
    it=1
    t=time()-ts
    status(it,re,t,cli,'start')
    while re>rt and it<max_it:
        Ux,Uy,p,Ux_,Uy_=_NS(Ux,Uy,p,gx,gy,nu,Nx,Ny,dx,dy,dt,bc,p_it)
        it+=1
        re=convergence(Ux,Uy,Ux_,Uy_)
        t=time()-ts
        status(it,re,t,cli,'run')
    F.update(Ux_,Uy_,p) # passing previous iteration to prevent errors
    F.stats(it,re,t)
    status(it,re,t,cli,'end')
    if cli: print('------------------------------------')

def convergence(Ux,Uy,Ux_,Uy_): # frobenius norm between velocity fields
    diff = norm(stack((Ux-Ux_,Uy-Uy_)).ravel())
    prev = norm(stack((Ux_,Uy_)).ravel())
    if prev < 1e-15: return inf
    return diff/prev