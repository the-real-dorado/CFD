from numba import njit,prange
from numpy import ones,zeros,linspace,meshgrid
from numpy import sum as arraysum
from matplotlib.pyplot import figure,contourf,axis,show
from time import time
from .utilities import importobj,in_cli,status

@njit(parallel=True) # numba go brrr
def _NS(vx,vy,p,gx,gy,rho,nu,Nx,Ny,dx,dy,dt): # solving the incompressible Navier-Stokes equations
    vx_,vy_=vx.copy(),vy.copy()               # with a poisson equation for the pressure field
    for _ in range(50):                       # using the intermediate velocity
        p_=p.copy()                           # by Jacobi iteration
        for i in prange(1,Nx-1):     
            for j in prange(1,Ny-1): 
                p[i,j]=( ( (p_[i+1,j]+p_[i-1,j])*dy**2 + (p_[i,j+1]+p_[i,j-1])*dx**2 ) / ( 2*(dx**2+dy**2) )
                        - ( ( rho*dx**2*dy**2 ) / ( 2*(dx**2+dy**2) ) )
                        * ( (1/dt)*((vx_[i+1,j]-vx_[i-1,j])/(2*dx)+(vy_[i,j+1]-vy_[i,j-1])/(2*dy))
                        - ((vx_[i+1,j]-vx_[i-1,j])/(2*dx))**2
                        - ((vy_[i,j+1]-vy_[i,j-1])/(2*dy))**2
                        -2*((vx_[i,j+1]-vx_[i,j-1])/(2*dy))*((vy_[i+1,j]-vy_[i-1,j])/(2*dx)) ) )
    for i in prange(1,Nx-1):     
        for j in prange(1,Ny-1):
            vx[i,j]=( vx_[i,j]
                   - vx_[i,j] * (dt/dx) * (vx_[i,j]-vx_[i-1,j]) - vy_[i,j] * (dt/dy) * (vx_[i,j]-vx_[i,j-1])
                   - (dt/(rho*2*dx)) * (p[i+1,j]-p[i-1,j])
                   + nu * (dt/dx**2) * (vx_[i+1,j]-2*vx_[i,j]+vx_[i-1,j]) + nu * (dt/dy**2) * (vx_[i,j+1]-2*vx_[i,j]+vx_[i,j-1])
                   + dt*gx )
            vy[i,j]=( vy_[i,j]
                   - vx_[i,j] * (dt/dx) * (vy_[i,j]-vy_[i-1,j]) - vy_[i,j] * (dt/dy) * (vy_[i,j]-vy_[i,j-1])
                   - (dt/(rho*2*dy)) * (p[i,j+1]-p[i,j-1])
                   + nu * (dt/dx**2) * (vy_[i+1,j]-2*vy_[i,j]+vy_[i-1,j]) + nu * (dt/dy**2) * (vy_[i,j+1]-2*vy_[i,j]+vy_[i,j-1])
                   + dt*gy )
    return vx,vy,p

def domain(file,dl,Nx,Ny,*,display=False): # domain for fluid
    if in_cli(): print('------------------------------------')
    print('Building domain')
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
    def __init__(self,domain,*,rho=1,nu=.1):
        self.domain,self.dl,self.Nx,self.Ny=domain
        self.rho=rho
        self.nu=nu
        self.x=linspace(0,(self.Nx-1)*self.dl,self.Nx)
        self.y=linspace(0,(self.Ny-1)*self.dl,self.Ny)
        self.X,self.Y=meshgrid(self.x,self.y)
    def update(self,vx,vy,p):
        self.vx=vx
        self.vy=vy
        self.p=p
    def stats(self,it,d,t):
        self.it=it
        self.d=d
        self.t=t

def simulation(fluid,bc,*,delta=1e-3,dt=1e-6,iterations=1000): # repeatedly calls calculation function
    F=fluid                                                    # until convergence
    Nx,Ny=F.Nx,F.Ny                                            # or max iterations
    dx=dy=F.dl
    vx,vy=zeros((Nx,Ny)),zeros((Nx,Ny))
    p=ones((Nx,Ny))
    rho,nu=F.rho,F.nu
    print('Initializing')
    vx,vy,p,gx,gy=bc(vx,vy,p,Nx,Ny)
    vx_,vy_=vx.copy(),vy.copy()
    vx,vy,p=_NS(vx,vy,p,gx,gy,rho,nu,Nx,Ny,dx,dy,dt)
    print(f'time-step: {dt}\nmax iterations: {iterations}\nmax delta-KE: {delta}')
    ts=time()
    cli=in_cli()
    if cli: print('------------------------------------')
    print(f'Solving')
    it=1
    d=convergence(vx,vy,vx_,vy_)
    t=time()-ts
    status(it,d,t,cli,'start')
    while d>delta and it<iterations:
        vx,vy,p,gx,gy=bc(vx,vy,p,Nx,Ny)
        vx_,vy_=vx.copy(),vy.copy()
        vx,vy,p=_NS(vx,vy,p,gx,gy,rho,nu,Nx,Ny,dx,dy,dt)
        it+=1
        d=convergence(vx,vy,vx_,vy_)
        t=time()-ts
        status(it,d,t,cli,'run')
    F.update(vx_,vy_,p) # passing previous iteration to prevent errors
    F.stats(it,d,t)
    status(it,d,t,cli,'end')
    if cli: print('------------------------------------')

def convergence(vx,vy,vx_,vy_): # fractional KE change
    d=0 # test if zero to prevent nan
    if arraysum(vx): d+=1-arraysum(vx_)**2/arraysum(vx)**2
    if arraysum(vy): d+=1-arraysum(vy_)**2/arraysum(vy)**2
    return abs(d)