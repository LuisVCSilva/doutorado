import numpy as np
import matplotlib.pyplot as plt

def passo(x, y, vx_meio, vy_meio, dt, beta, k, F, m):
    fx = np.random.normal(0, 1.0)
    fy = np.random.normal(0, 1.0)
    x_prox = x + dt * vx_meio
    y_prox = y + dt * vy_meio
    ax = (-beta * vx_meio - k * x_prox + F * fx) / m
    ay = (-beta * vy_meio - k * y_prox + F * fy) / m
    vx_meio_prox = vx_meio + dt * ax
    vy_meio_prox = vy_meio + dt * ay
    return x_prox, y_prox, vx_meio_prox, vy_meio_prox

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
    x[i+1], y[i+1], vx_meio[i+1], vy_meio[i+1] = passo(
        x[i], y[i], vx_meio[i], vy_meio[i], dt, beta, k, F, m
    )

plt.figure(figsize=(6,6))
plt.plot(x, y, lw=0.5)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Trajetória da Partícula - Leapfrog com Função passo')
plt.axis('equal')
plt.grid(True)
plt.show()

