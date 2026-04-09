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
Nx = 21
Nt = 100

dx = L / (Nx - 1)
dt = T / Nt

x = np.linspace(0, L, Nx)
t = np.linspace(0, T, Nt)

# ---------------------------------
# D1Q3
# ---------------------------------
c = np.array([0, 1, -1])
w = np.array([2/3, 1/6, 1/6])

# velocidade de lattice
cs = dx / dt

# relaxação
tau = 0.51

# ---------------------------------
# Distribuições
# ---------------------------------
f = np.zeros((3, Nx))
f_new = np.zeros((3, Nx))

# ---------------------------------
# Condição inicial
# ---------------------------------
p = np.cos(x) * np.sin(2*np.pi*x)
u = np.cos(2*np.pi*x)

# ---------------------------------
# Inicialização em equilíbrio
# ---------------------------------
for i in range(3):
    f[i,:] = w[i] * (p + alpha*c[i]*u/(cs**2))

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
    u = np.sum(c[:,None] * f, axis=0)

    # equilíbrio
    feq = np.zeros_like(f)
    for i in range(3):
        feq[i,:] = w[i] * (p + alpha*c[i]*u/(cs**2))

    # colisão
    f = f - (1/tau)*(f-feq)

    # streaming
    f_new[0,:] = f[0,:]
    f_new[1,1:] = f[1,:-1]
    f_new[2,:-1] = f[2,1:]

    # contorno pressão
    p_left = 0.5
    p_right = 1.8

    f_new[:,0] = w * p_left
    f_new[:,-1] = w * p_right

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

plt.title('LBM D1Q3 para Equação da Onda')
plt.show()
