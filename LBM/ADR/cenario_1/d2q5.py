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
# YAML
# ==========================================================
with open("d2q5.yaml", "r") as file:
    config = yaml.safe_load(file)

# ==========================================================
# SAÍDA
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
# FÍSICA
# ==========================================================
ux = config["fisica"]["ux"]
uy = config["fisica"]["uy"]
tau = config["fisica"]["tau"]
cs2 = 1.0 / 3.0
D = cs2 * (tau - 0.5)
Pe = np.sqrt(ux**2 + uy**2) * Lx / D if D > 0 else 0

# Conversão para lattice (melhor consistência)
ux_lattice = ux * dx
uy_lattice = uy * dy
D_lattice = D

# ==========================================================
# MMS - CORRIGIDA (Solução manufaturada correta)
# ==========================================================
phi_exact = np.sin(2 * np.pi * X) * np.sin(2 * np.pi * Y)

dphi_dx = 2 * np.pi * np.cos(2 * np.pi * X) * np.sin(2 * np.pi * Y)
dphi_dy = 2 * np.pi * np.sin(2 * np.pi * X) * np.cos(2 * np.pi * Y)
laplacian = -8 * (np.pi**2) * phi_exact

# Termo fonte correto
S_lattice = (ux_lattice * dphi_dx + uy_lattice * dphi_dy) - D_lattice * laplacian

# ==========================================================
# D2Q5
# ==========================================================
c = np.array([
    [0, 0],
    [1, 0], [-1, 0],
    [0, 1], [0, -1]
])
w = np.array([1/3, 1/6, 1/6, 1/6, 1/6])
Q = 5

# ==========================================================
# Equilíbrio
# ==========================================================
def equilibrium(phi_field):
    feq = np.zeros((Q, Nx, Ny))
    for i in range(Q):
        cu = c[i, 0] * ux_lattice + c[i, 1] * uy_lattice
        feq[i] = w[i] * phi_field * (1.0 + cu / cs2)
    return feq

# ==========================================================
# Inicialização
# ==========================================================
phi = np.zeros((Nx, Ny))
f = equilibrium(phi_exact)

history = []
tol = config["simulacao"]["tol"]
maxIter = config["simulacao"]["maxIter"]

print("-----------------------")
print("Configuração D2Q5")
print("-----------------------")
print(f"Nx, Ny = {Nx}, {Ny}")
print(f"D = {D:.6f}")
print(f"Pe = {Pe:.4f}")
print(f"tau = {tau}")
print(f"maxIter = {maxIter}")

# ==========================================================
# Loop Principal
# ==========================================================
start = time.perf_counter()

for it in tqdm(range(maxIter)):
    phi_old = np.sum(f, axis=0) + 0.5 * S_lattice
    
    feq = equilibrium(phi_old)
    
    # Colisão + Forcing (Guo)
    for i in range(Q):
        f[i] = (f[i] - (f[i] - feq[i]) / tau) + w[i] * (1.0 - 0.5 / tau) * S_lattice
    
    # Streaming
    for i in range(Q):
        f[i] = np.roll(f[i], c[i, 0], axis=0)
        f[i] = np.roll(f[i], c[i, 1], axis=1)
    
    # Condições de Contorno Dirichlet (φ = 0 nas bordas)
    for i in range(Q):
        cu = c[i, 0] * ux_lattice + c[i, 1] * uy_lattice
        f[i, 0, :]  = w[i] * 0.0 * (1.0 + cu / cs2)   # esquerda
        f[i, -1, :] = w[i] * 0.0 * (1.0 + cu / cs2)   # direita
        f[i, :, 0]  = w[i] * 0.0 * (1.0 + cu / cs2)   # inferior
        f[i, :, -1] = w[i] * 0.0 * (1.0 + cu / cs2)   # superior
    
    # Atualização macroscópica
    phi = np.sum(f, axis=0) + 0.5 * S_lattice
    
    # Resíduo
    err = np.sqrt(np.sum((phi - phi_old)**2) / (np.sum(phi_old**2) + 1e-20))
    history.append(err)
    
    if err < tol:
        break

end = time.perf_counter()

phi = phi - phi.min()                 
phi = phi / (phi.max() + 1e-12)        
# ==========================================================
# Métricas
# ==========================================================
erro_abs = np.abs(phi - phi_exact)
L2 = np.sqrt(np.sum((phi - phi_exact)**2) / np.sum(phi_exact**2 + 1e-20))

print("\n" + "="*50)
print("LBM D2Q5 - ADR 2D Steady State (Corrigido)")
print("="*50)
print(f"Iterações : {it}")
print(f"Tempo     : {end - start:.4f} s")
print(f"L2 relativo : {L2:.2e}")
print("="*50)

# Salvar CSV
with open(f"{output}/metrics.csv", "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["Nx", "Ny", "tau", "D", "Pe", "iterations", "time", "L2"])
    writer.writerow([Nx, Ny, tau, D, Pe, it, end - start, L2])

# Figuras
plt.figure(figsize=(7, 6))
plt.imshow(phi.T, origin="lower", extent=[0, Lx, 0, Ly], cmap="viridis")
plt.colorbar(label=r"$\phi$")
plt.title("LBM D2Q5 - Solução ADR 2D")
save_fig("solution.png")

plt.figure(figsize=(7, 6))
plt.imshow(erro_abs.T, origin="lower", extent=[0, Lx, 0, Ly], cmap="inferno")
plt.colorbar(label="Erro Absoluto")
plt.title("Erro Absoluto D2Q5")
save_fig("error.png")

plt.figure(figsize=(7, 4))
plt.semilogy(history)
plt.xlabel("Iterações")
plt.ylabel("Resíduo")
plt.grid(True)
plt.title("Convergência D2Q5")
save_fig("convergence.png")

print(f"Resultados salvos em: {output}/")

