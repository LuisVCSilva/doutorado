import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------
# Parâmetros físicos
# ---------------------------------
L = 1.0
T = 1.0
alpha = np.sqrt(2)

# ---------------------------------
# Malha
# ---------------------------------
Nx = 101
Nt = 500

dx = L / (Nx - 1)
dt = T / Nt

x = np.linspace(0, L, Nx)
t = np.linspace(0, T, Nt)

# ---------------------------------
# D1Q5
# ---------------------------------
c = np.array([0, 1, -1, 2, -2])

# pesos corretos
w = np.array([0.5, 1/6, 1/6, 1/12, 1/12])

# velocidade lattice
cs = dx / dt

# relaxação mínima dissipativa
tau = 0.50001

# ---------------------------------
# Distribuições
# ---------------------------------
f = np.zeros((5, Nx))
f_new = np.zeros((5, Nx))

# ---------------------------------
# Condição inicial
# ---------------------------------
p = np.cos(x) * np.sin(2*np.pi*x)
u = np.cos(2*np.pi*x)

# ---------------------------------
# Inicialização correta
# ---------------------------------
for i in range(5):
    f[i,:] = w[i] * (
        p
        + c[i]*u/(cs**2)
        + ((c[i]**2 - cs**2)*p)/(2*cs**4)
    )

# ---------------------------------
# Histórico
# ---------------------------------
history = np.zeros((Nt, Nx))

# ---------------------------------
# Evolução temporal
# ---------------------------------
for n in range(Nt):

    # macroscópicas
    p = np.sum(f, axis=0)
    u = cs * np.sum(c[:,None] * f, axis=0)

    # equilíbrio
    feq = np.zeros_like(f)
    for i in range(5):
        feq[i,:] = w[i] * (
            p
            + c[i]*u/(cs**2)
            + ((c[i]**2 - cs**2)*p)/(2*cs**4)
        )

    # colisão
    f = f - (1/tau)*(f-feq)

    # ---------------------------------
    # streaming
    # ---------------------------------
    f_new[:] = 0

    # c=0
    f_new[0,:] = f[0,:]

    # c=+1
    f_new[1,1:] = f[1,:-1]

    # c=-1
    f_new[2,:-1] = f[2,1:]

    # c=+2
    f_new[3,2:] = f[3,:-2]

    # c=-2
    f_new[4,:-2] = f[4,2:]

    # ---------------------------------
    # Condições de contorno
    # ---------------------------------
    p_left = 0.5
    p_right = 1.8

    for i in range(5):
        f_new[i,0] = w[i] * p_left
        f_new[i,-1] = w[i] * p_right

    f = f_new.copy()

    history[n,:] = np.sum(f, axis=0)

# ---------------------------------
# Valor pedido
# ---------------------------------
ix = int(0.5/dx)
print("p(0.5,1) =", history[-1, ix])

# ---------------------------------
# Plot 3D
# ---------------------------------
X, Y = np.meshgrid(x, t)

fig = plt.figure(figsize=(10,7))
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(X, Y, history)

ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('p(x,t)')

plt.title('LBM D1Q5 Corrigido - Equação da Onda')

plt.show()
