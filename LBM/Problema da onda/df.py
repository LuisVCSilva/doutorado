import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# -----------------------------
# Parâmetros físicos
# -----------------------------
L = 1.0
T = 1.0
alpha = np.sqrt(2)

# -----------------------------
# Parâmetros numéricos
# -----------------------------
h = 0.05
k = 0.01

m = int(L / h)
N = int(T / k)

x = np.linspace(0, L, m + 1)
t = np.linspace(0, T, N + 1)

lam = alpha * k / h

print(f"lambda = {lam:.4f}")

if lam > 1:
    raise ValueError("Método instável: lambda > 1")

# -----------------------------
# Inicialização
# -----------------------------
w = np.zeros((m + 1, N + 1))

# -----------------------------
# Condições de contorno
# -----------------------------
w[0, :] = 0.5
w[m, :] = 1.8

# -----------------------------
# Condição inicial
# -----------------------------
f = lambda x: np.cos(x) * np.sin(2 * np.pi * x)
g = lambda x: np.cos(2 * np.pi * x)

w[:, 0] = f(x)
w[0, 0] = 0.5
w[m, 0] = 1.8

# -----------------------------
# Primeira camada temporal
# -----------------------------
for i in range(1, m):
    w[i, 1] = (
        (1 - lam**2) * f(x[i])
        + (lam**2 / 2) * (f(x[i + 1]) + f(x[i - 1]))
        + k * g(x[i])
    )

# -----------------------------
# Evolução temporal
# -----------------------------
for j in range(1, N):
    for i in range(1, m):
        w[i, j + 1] = (
            2 * (1 - lam**2) * w[i, j]
            + lam**2 * (w[i + 1, j] + w[i - 1, j])
            - w[i, j - 1]
        )

# -----------------------------
# Valor pedido
# -----------------------------
ix = int(0.5 / h)
jt = int(1.0 / k)

print(f"p(0.5,1) = {w[ix,jt]:.6f}")

# -----------------------------
# Gráfico 3D
# -----------------------------
X, Tm = np.meshgrid(t, x)

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(Tm, X, w, edgecolor='none')

ax.set_xlabel('Tempo t')
ax.set_ylabel('Posição x')
ax.set_zlabel('Pressão p(x,t)')

plt.title('Propagação da pressão na mangueira')
plt.show()
