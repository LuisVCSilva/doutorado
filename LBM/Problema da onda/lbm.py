import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------
# Parâmetros físicos
# ---------------------------------
L = 1.0
T = 1.0

Nx = 21
Nt = 100

dx = L/(Nx-1)
dt = T/Nt

x = np.linspace(0, L, Nx)
t = np.linspace(0, T, Nt)

# ---------------------------------
# D1Q3
# ---------------------------------
c = np.array([0, 1, -1])
w = np.array([2/3, 1/6, 1/6])

tau = 1.0

# ---------------------------------
# Distribuições
# ---------------------------------
f = np.zeros((3, Nx))
f_new = np.zeros((3, Nx))

# ---------------------------------
# Condição inicial
# ---------------------------------
p = np.cos(x) * np.sin(2*np.pi*x)

for i in range(3):
    f[i,:] = w[i]*p

# ---------------------------------
# Histórico
# ---------------------------------
history = np.zeros((Nt, Nx))

# ---------------------------------
# Evolução temporal
# ---------------------------------
for n in range(Nt):

    # colisão
    p = np.sum(f, axis=0)

    feq = np.zeros_like(f)
    for i in range(3):
        feq[i,:] = w[i]*p

    f = f - (1/tau)*(f-feq)

    # streaming
    f_new[0,:] = f[0,:]
    f_new[1,1:] = f[1,:-1]
    f_new[2,:-1] = f[2,1:]

    # contorno
    f_new[:,0] = w*0.5
    f_new[:,-1] = w*1.8

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

plt.title('LBM D1Q3 - Propagação de pressão')
plt.show()
