import numpy as np
import matplotlib.pyplot as plt

T = 100         
dt = 1e-3       
N = int(T / dt) 
t = np.linspace(0, T, N)

beta = 1.0 
F = 10.0    
k = 0.5     
m = 1.0    

x = np.zeros(N)
y = np.zeros(N)
vx_meio = np.zeros(N)  
vy_meio = np.zeros(N)

x[0] = 1.0
y[0] = 0.0
vx0 = 0.0
vy0 = 1.0


fx0 = np.random.normal(0, 1.0)
fy0 = np.random.normal(0, 1.0)
ax0 = (-beta * vx0 - k * x[0] + F * fx0) / m
ay0 = (-beta * vy0 - k * y[0] + F * fy0) / m

vx_meio[0] = vx0 + 0.5 * dt * ax0
vy_meio[0] = vy0 + 0.5 * dt * ay0

for i in range(N - 1):

    x[i+1] = x[i] + dt * vx_meio[i]
    y[i+1] = y[i] + dt * vy_meio[i]


    fx = np.random.normal(0, 1.0)
    fy = np.random.normal(0, 1.0)


    ax = (-beta * vx_meio[i] - k * x[i+1] + F * fx) / m
    ay = (-beta * vy_meio[i] - k * y[i+1] + F * fy) / m


    vx_meio[i+1] = vx_meio[i] + dt * ax
    vy_meio[i+1] = vy_meio[i] + dt * ay

plt.figure(figsize=(6,6))
plt.plot(x, y, lw=0.5)
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')
plt.grid(True)
plt.show()
