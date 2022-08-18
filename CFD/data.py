from matplotlib.pyplot import figure,contourf,streamplot,colorbar,axis,title,savefig
from matplotlib.colors import colorConverter,LinearSegmentedColormap
from .utilities import in_cli

# custom colormap for overlaying domain
color_w=colorConverter.to_rgba('white',alpha=0)
color_b=colorConverter.to_rgba('black',alpha=1)
cmap_alpha=LinearSegmentedColormap.from_list('cmap_alpha',[color_w,color_b],2)

def plot(fluid,plots=['velocity','streamlines'],*,path='',info=False):
   name='plot'
   for _ in plots: name=name+'-'+_
   if in_cli(): print(f'Rendering: {name}')
   F=fluid
   figure(figsize=((F.Nx-1)/(F.Ny-1)*10,10),tight_layout=True,facecolor='w')
   if 'pressure' in plots: 
      contourf(F.X,F.Y,F.p.T,cmap='bwr',levels=20)
      if info: title('\nPressure [Pa]\n',fontdict={'fontsize':30})
   if 'velocity' in plots: 
      contourf(F.X,F.Y,(F.vx.T**2+F.vy.T**2)**0.5,cmap='Blues',levels=20)
      if info: title('\nVelocity [m/s]\n',fontdict={'fontsize':30})
   if info: colorbar()
   if 'streamlines' in plots: streamplot(F.X[2:-3][2:-3],F.Y[2:-3][2:-3],F.vx.T[2:-3][2:-3],F.vy.T[2:-3][2:-3],color='k',density=2,arrowstyle='->')
   contourf(F.x,F.y,-F.domain.T,cmap=cmap_alpha)
   axis('scaled')
   axis('off')
   savefig(f'{path}{name}_it-{F.it}_d-{F.d:.5f}_t-{F.t//60:.0f}m{F.t%60:.0f}s.png',dpi=200,bbox_inches='tight')
   if in_cli(): print('------------------------------------')

def export(fluid,path=''):
   name='export'
   F=fluid
   print(f'Exporting')
   with open(f'{path}{name}_it-{F.it}_d-{F.d:.5f}_t-{F.t//60:.0f}m{F.t%60:.0f}s.csv','w+') as f:
      f.write('x,y,vx,vy,p\n')
      for i in range(F.Nx):
         for j in range(F.Ny):
            f.write(f'{i*F.dl},{j*F.dl},{F.vx[i,j]},{F.vy[i,j]},{F.p[i,j]}\n')
   if in_cli(): print('------------------------------------')