# simulation requires a function
# which enforces boundary conditions 
# on each iteration

def flow(domain,gx,gy): # flow with obstructions
    def bc(vx,vy,p,Nx,Ny):
        for i in range(Nx):
            for j in range(Ny):
                if not domain[0][i,j]: # obstructions
                    p[i+1,j+1]=p[i-1,j-1]=p[i-1,j+1]=p[i+1,j-1]
                    vx[i,j]=vy[i,j]=0
        vx[:,0]=vx[:,1]
        vy[:,0]=vy[:,1]
        p[:,0]=p[:,1]
        vx[:,-1]=vx[:,-2]
        vy[:,-1]=vy[:,-2]
        p[:,-1]=p[:,-2]
        vx[0,:]=vx[1,:]
        vy[0,:]=vy[1,:]
        p[0,:]=p[1,:]
        vx[-1,:]=vx[-2,:]
        vy[-1,:]=vy[-2,:]
        p[-1,:]=p[-2,:]
        return vx,vy,p,gx,gy
    return bc

def lid(vx): # standard lid-driven cavity problem
    _=vx
    def bc(vx,vy,p,Nx,Ny):
        p[-1,:]=p[-2,:]
        p[0,:]=p[1,:]
        p[:,0]=p[:,1]
        p[:,-1]=0
        vx[-1,:]=vy[-1,:]=0
        vx[0,:]=vy[0,:]=0
        vx[:,0]=vy[:,0]=0
        vx[:,-1]=_
        vy[:,-1]=0
        return vx,vy,p,0,0
    return bc

def channel(gx): # channel flow
    def bc(vx,vy,p,Nx,Ny):
        vx[:,0]=0
        vy[:,0]=0
        p[:,0]=p[:,1]=p[:,2]
        vx[:,-1]=0
        vy[:,-1]=0
        p[:,-1]=p[:,-2]=p[:,-3]
        vx[0,:]=vx[1,:]
        vy[0,:]=vy[1,:]
        p[0,:]=p[1,:]
        vx[-1,:]=vx[-2,:]
        vy[-1,:]=vy[-2,:]
        p[-1,:]=p[-2,:]
        return vx,vy,p,gx,0
    return bc