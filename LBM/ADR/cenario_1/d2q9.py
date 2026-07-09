import numpy as np
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import csv
import time
import os
from tqdm import tqdm
import yaml


# ==========================================================
# LEITURA DO YAML
# ==========================================================

with open("d2q9.yaml", "r") as file:
    config = yaml.safe_load(file)


# ==========================================================
# FUNÇÃO PARA SALVAR FIGURAS
# ==========================================================

output = config["saida"]["diretorio"]
os.makedirs(output, exist_ok=True)


def save_fig(filename):
    plt.tight_layout()
    plt.savefig(f"{output}/{filename}", dpi=300, bbox_inches="tight")
    plt.close()


# ==========================================================
# DOMÍNIO
# ==========================================================

Nx = config["dominio"]["Nx"]
Ny = config["dominio"]["Ny"]

Lx = config["dominio"]["Lx"]
Ly = config["dominio"]["Ly"]

dx = Lx / (Nx - 1)
dy = Ly / (Ny - 1)

x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)

X, Y = np.meshgrid(x, y, indexing="ij")


# ==========================================================
# PARÂMETROS FÍSICOS
# ==========================================================

ux = config["fisica"]["ux"]
uy = config["fisica"]["uy"]
tau = config["fisica"]["tau"]


cs2 = 1.0 / 3.0

# Difusão física do LBM
D = cs2 * (tau - 0.5)

# Número de Peclet no contínuo
Pe = np.sqrt(ux**2 + uy**2) * Lx / D


# ==========================================================
# PARÂMETROS DO LATTICE
# ==========================================================

# Conversão para unidades lattice
ux_lattice = ux * dx
uy_lattice = uy * dy

D_lattice = cs2 * (tau - 0.5)


# ==========================================================
# SOLUÇÃO MANUFATURADA (MMS)
# ==========================================================

kx = np.pi / (Nx - 1)
ky = np.pi / (Ny - 1)


X_lat, Y_lat = np.meshgrid(
    np.arange(Nx),
    np.arange(Ny),
    indexing="ij"
)


phi_exact = (
    np.sin(kx * X_lat) *
    np.sin(ky * Y_lat)
)


# Derivadas analíticas

dphi_dx = (
    kx *
    np.cos(kx * X_lat) *
    np.sin(ky * Y_lat)
)


dphi_dy = (
    ky *
    np.sin(kx * X_lat) *
    np.cos(ky * Y_lat)
)


laplacian = (
    -(kx**2 + ky**2)
    * phi_exact
)


# Fonte ADR:
# u.grad(phi) = D laplacian(phi) + S
# Portanto:
# S = u.grad(phi) - D laplacian(phi)

S_lattice = (
    ux_lattice * dphi_dx
    +
    uy_lattice * dphi_dy
    -
    D_lattice * laplacian
)


# ==========================================================
# PARÂMETROS DA SIMULAÇÃO
# ==========================================================

tol = config["simulacao"]["tol"]
maxIter = config["simulacao"]["maxIter"]


print("Configuração carregada:")
print("-----------------------")
print(f"Nx, Ny       = {Nx}, {Ny}")
print(f"Lx, Ly       = {Lx}, {Ly}")
print(f"ux, uy       = {ux}, {uy}")
print(f"tau          = {tau}")
print(f"D            = {D}")
print(f"Pe           = {Pe}")
print(f"tol          = {tol}")
print(f"maxIter      = {maxIter}")
print(f"Saída        = {output}")

# ==========================================================
c = np.array([
    [0, 0], [1, 0], [0, 1], [-1, 0], [0, -1],
    [1, 1], [-1, 1], [-1, -1], [1, -1]
])
w = np.array([4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36])
Q = 9

# Inicialização
phi = phi_exact.copy()
f = np.zeros((Q, Nx, Ny))

def get_equilibrium(phi_field):
    feq = np.zeros((Q, Nx, Ny))
    for i in range(Q):
        cu = c[i, 0] * ux_lattice + c[i, 1] * uy_lattice
        feq[i] = w[i] * phi_field * (1.0 + cu / cs2)
    return feq

f = get_equilibrium(phi)

# ==========================================================
tol = 1e-12
maxIter = 50000
history = []
start = time.perf_counter()

for it in tqdm(range(maxIter)):
    phi_old = np.sum(f, axis=0) + 0.5 * S_lattice
    
    # Colisão BGK + Fonte de Guo
    feq = get_equilibrium(phi_old)
    for i in range(Q):
        f[i] = f[i] - (f[i] - feq[i]) / tau + w[i] * (1.0 - 0.5 / tau) * S_lattice
        
    # Streaming
    for i in range(Q):
        f[i] = np.roll(f[i], c[i, 0], axis=0)
        f[i] = np.roll(f[i], c[i, 1], axis=1)
        
    # Condições de Contorno Dirichlet (Injeção de equilíbrio)
    for i in range(Q):
        cu = c[i, 0] * ux_lattice + c[i, 1] * uy_lattice
        f[i, 0, :]  = w[i] * phi_exact[0, :]  * (1.0 + cu / cs2)
        f[i, -1, :] = w[i] * phi_exact[-1, :] * (1.0 + cu / cs2)
        f[i, :, 0]  = w[i] * phi_exact[:, 0]  * (1.0 + cu / cs2)
        f[i, :, -1] = w[i] * phi_exact[:, -1] * (1.0 + cu / cs2)

    # Atualização macroscópica com correção de Guo
    phi = np.sum(f, axis=0) + 0.5 * S_lattice
    
    # Resíduo para histórico de convergência
    err = np.sqrt(np.sum((phi - phi_old)**2) / (np.sum(phi_old**2) + 1e-15))
    history.append(err)
    if err < tol:
        break

end = time.perf_counter()

# ==========================================================
erro = np.abs(phi - phi_exact)
L2 = np.sqrt(np.sum((phi - phi_exact)**2) / np.sum(phi_exact**2))

print("\n" + "="*50)
print("LBM D2Q9 - ADR 2D Steady State")
print("="*50)
print(f"Iterações : {it}")
print(f"Tempo     : {end - start:.4f} s")
print(f"D         : {D:.6f}")
print(f"Pe        : {Pe:.4f}")
print(f"L2 relativo: {L2:.2e}")
print("="*50)

# Salvar resultados
with open(output + "/metrics.csv", "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["Nx", "Ny", "tau", "D", "Pe", "iterations", "time", "L2"])
    writer.writerow([Nx, Ny, tau, D, Pe, it, end - start, L2])

with open(output + "/solution.csv", "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["x", "y", "phi_LBM", "phi_exact", "error"])
    for i in range(Nx):
        for j in range(Ny):
            writer.writerow([X[i, j], Y[i, j], phi[i, j], phi_exact[i, j], erro[i, j]])

# Figuras
plt.figure(figsize=(7, 6))
plt.imshow(phi.T, origin="lower", extent=[0, Lx, 0, Ly], cmap="viridis")
plt.colorbar(label=r"$\phi$")
plt.title("LBM D2Q9 - Solução ADR 2D")
save_fig("solution.png")

plt.figure(figsize=(7, 6))
plt.imshow(erro.T, origin="lower", extent=[0, Lx, 0, Ly], cmap="inferno")
plt.colorbar(label="Erro Absoluto")
plt.title("Erro Absoluto")
save_fig("error.png")

plt.figure(figsize=(7, 4))
plt.semilogy(history)
plt.xlabel("Iterações")
plt.ylabel("Resíduo")
plt.grid(True)
plt.title("Convergência D2Q9")
save_fig("convergence.png")

print(f"Resultados salvos em: {output}/")
