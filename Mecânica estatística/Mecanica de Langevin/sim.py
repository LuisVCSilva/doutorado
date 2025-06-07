import numpy as np
import matplotlib.pyplot as plt

# Parâmetros
T = 100
dt = 1e-3
N = int(T / dt)
t = np.linspace(0, T, N)

# Parâmetros físicos
beta = 1.0      # Coef. de dissipação
F = 10.0        # Intensidade do ruído
k = 0.5         # Constante da mola (Hooke)
N_exps = 100    # Número de simulações para média
m = 1.0         # Massa da partícula


vx_med = np.zeros(N)
vy_med = np.zeros(N)
v2_med = np.zeros(N)
energia_med = np.zeros(N)

for _ in range(N_exps):

    x = np.zeros(N)
    y = np.zeros(N)
    vx = np.zeros(N)
    vy = np.zeros(N)


    x[0] = 1.0
    y[0] = 0.0
    vx[0] = 0.0
    vy[0] = 0.0

    for i in range(N - 1):

        fx = np.random.normal(0, 1.0)
        fy = np.random.normal(0, 1.0)


        vx[i+1] = vx[i] + dt * (-beta * vx[i] - k * x[i] + F * fx)
        vy[i+1] = vy[i] + dt * (-beta * vy[i] - k * y[i] + F * fy)


        x[i+1] = x[i] + dt * vx[i]
        y[i+1] = y[i] + dt * vy[i]


    vx_med += vx
    vy_med += vy
    v2_med += vx**2 + vy**2
    energia_med += 0.5 * m * (vx**2 + vy**2)

vx_med /= N_exps
vy_med /= N_exps
v2_med /= N_exps
energia_med /= N_exps

plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.plot(t, vx_med, label='⟨vₓ⟩')
plt.plot(t, vy_med, label='⟨vᵧ⟩')
plt.xlabel('Tempo')
plt.ylabel('Velocidade média')
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(t, energia_med, label='⟨E⟩ = energia cinética média')
plt.xlabel('Tempo')
plt.ylabel('Energia')
plt.legend()
plt.tight_layout()
plt.show()
