import matplotlib.pyplot as plt
from scipy import *
from scipy import integrate
from scipy.integrate import ode, odeint
import numpy as np

plt.close()
fig,ax = plt.subplots(1,1,figsize=(8,8))
ax.set_facecolor(plt.cm.gray(.95))
#ax.axis('equal')
## Vector field function
def vf1(arg,t, a,b):
  x = arg[0]
  y = arg[1]
  dx=np.zeros(2)
  dx[0]=-x +a*y + x**2*y
  dx[1]=b - a*y - x**2*y
  return dx
#Integrate under curve
# Initial condition
xx0 = np.array([.25,1.5])
a = .08
b = .6
# Solve
t = np.linspace(0, 15, 1000)
ys = odeint(vf1, xx0, t, args=(a,b))

#Vector field
x1,y1 = np.linspace(0,3,200),np.linspace(0,b/a,200)
X,Y = np.meshgrid(x1,y1 )
U = -X +a*Y + X**2*Y
V = b - a*Y - X**2*Y

start = [[0,2],[.8,0]]

strm = ax.streamplot( X,Y,U, V,linewidth=.2)
strmS = ax.streamplot(x1,y1, U, V, start_points=start, color="crimson", linewidth=1)

ax.plot(ys[:,0],ys[:,1],label=r'para $x_0 = ({0}, {1})$'.format(xx0[0],xx0[1]),lw=.5)
ax.plot(b,b/(a+b**2), 'k*', label='ponto fixo em $({0}, {1})$'.format(b,np.round(b/(a+b**2),2)))
ax.plot(x1,b/(a+x1**2), 'k-.', lw=.8, label=r'y-nullcline $ \frac{b}{(a + x^2)}$')
ax.plot(x1,x1/(a+x1**2), 'k:',lw=1, label=r'x-nullcline $ \frac{x}{(a + x^2)}$')
#plt.xlim([-1,3])
plt.ylim([0,b/a+.2])
plt.xlabel(r"$ x$",size=14)
plt.ylabel(r"$ y$",size=14)
plt.legend()
plt.grid(True)
plt.show()
