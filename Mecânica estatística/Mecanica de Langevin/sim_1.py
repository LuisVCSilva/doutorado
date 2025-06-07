import numpy as np
import matplotlib.pyplot as plt


T = 100          # Tempo total
dt = 1e-3        # Passo de tempo
N = int(T / dt)  # Número de passos
t = np.linspace(0, T, N)


beta = 1.0   # Coeficiente de dissipação
F = 10.0     # Intensidade do ruído
k = 0.5      # Constante da mola
m = 1.0      # Massa


x = np.zeros(N)
y = np.zeros(N)
vx = np.zeros(N)
vy = np.zeros(N)


x[0] = 1.0
y[0] = 0.0
vx[0] = 0.0
vy[0] = 1.0


for i in range(N - 1):
    fx = np.random.normal(0, 1.0)
    fy = np.random.normal(0, 1.0)

    vx[i+1] = vx[i] + dt * (-beta * vx[i] - k * x[i] + F * fx)
    vy[i+1] = vy[i] + dt * (-beta * vy[i] - k * y[i] + F * fy)


    x[i+1] = x[i] + dt * vx[i]
    y[i+1] = y[i] + dt * vy[i]

plt.figure(figsize=(6,6))
plt.plot(x, y, lw=0.5)
plt.title('Trajetória da Partícula no Plano (x, y)')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')
plt.grid(True)
plt.show()

