import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import time
from tqdm import tqdm

# ==========================================================
# Configuração
# ==========================================================
output = "dados/d2q9_adr_2d_results"
os.makedirs(output, exist_ok=True)

# Carregue seu config ou defina aqui
Nx = 100
Ny = 100
Lx = 1.0
Ly = 1.0

ux = 0.0      # mude conforme necessário
uy = 0.0
tau = 0.8     # ajuste conforme o D desejado

# ==========================================================
# Parâmetros da malha
# ==========================================================
dx = Lx / (Nx - 1)
dy = Ly / (Ny - 1)
x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)
X, Y = np.meshgrid(x, y, indexing="ij")

# ==========================================================
# Parâmetros físicos
# ==========================================================
cs2 = 1.0 / 3.0
D = cs2 * (tau - 0.5)                    # Difusividade física
Pe = np.sqrt(ux**2 + uy**2) * Lx / D if D > 0 else 0

# Conversão para lattice
ux_lattice = ux * dx
uy_lattice = uy * dy
D_lattice = cs2 * (tau - 0.5)

print(f"D = {D:.6f} | Pe = {Pe:.4f} | tau = {tau}")

# ==========================================================
# Solução Manufaturada Correta (MMS)
# ==========================================================
phi_exact = np.sin(2 * np.pi * X) * np.sin(2 * np.pi * Y)

dphi_dx = 2 * np.pi * np.cos(2 * np.pi * X) * np.sin(2 * np.pi * Y)
dphi_dy = 2 * np.pi * np.sin(2 * np.pi * X) * np.cos(2 * np.pi * Y)
laplacian = -8 * (np.pi**2) * phi_exact

# Termo fonte S (para que phi_exact seja solução exata)
S_lattice = (ux_lattice * dphi_dx + uy_lattice * dphi_dy) - D_lattice * laplacian

# ==========================================================
# Lattice D2Q9
# ==========================================================
c = np.array([
    [0, 0], [1, 0], [0, 1], [-1, 0], [0, -1],
    [1, 1], [-1, 1], [-1, -1], [1, -1]
])
w = np.array([4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36])
Q = 9

# Inicialização
phi = np.zeros((Nx, Ny))
f = np.zeros((Q, Nx, Ny))

def get_equilibrium(phi_field, ux_l, uy_l):
    feq = np.zeros((Q, Nx, Ny))
    for i in range(Q):
        cu = c[i, 0] * ux_l + c[i, 1] * uy_l
        feq[i] = w[i] * phi_field * (1.0 + cu / cs2)
    return feq

# Inicializa com equilíbrio
f = get_equilibrium(phi_exact, ux_lattice, uy_lattice)

# ==========================================================
# Simulação
# ==========================================================
tol = 1e-10
maxIter = 50000
history = []

start = time.perf_counter()

for it in tqdm(range(maxIter)):
    phi_old = np.sum(f, axis=0) + 0.5 * S_lattice
    
    # Colisão + Forcing (Guo)
    feq = get_equilibrium(phi_old, ux_lattice, uy_lattice)
    for i in range(Q):
        f[i] = f[i] - (f[i] - feq[i]) / tau + w[i] * (1.0 - 0.5 / tau) * S_lattice
    
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
    
    # Residuo
    err = np.sqrt(np.sum((phi - phi_old)**2) / (np.sum(phi_old**2) + 1e-20))
    history.append(err)
    
    if err < tol:
        break
phi = phi - phi.min()                  
phi = phi / (phi.max() + 1e-12)        
end = time.perf_counter()

# ==========================================================
# Métricas e Saída
# ==========================================================
erro_abs = np.abs(phi - phi_exact)
L2 = np.sqrt(np.sum((phi - phi_exact)**2) / np.sum(phi_exact**2))

print("\n" + "="*60)
print("LBM D2Q9 - ADR 2D Steady State (Corrigido)")
print("="*60)
print(f"Iterações : {it}")
print(f"Tempo     : {end - start:.4f} s")
print(f"L2 relativo : {L2:.2e}")
print("="*60)

# Salvar resultados
with open(f"{output}/metrics.csv", "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["Nx", "Ny", "tau", "D", "Pe", "iterations", "time", "L2"])
    writer.writerow([Nx, Ny, tau, D, Pe, it, end-start, L2])

# Plot Solução
plt.figure(figsize=(7, 6))
plt.imshow(phi.T, origin="lower", extent=[0, Lx, 0, Ly], cmap="viridis")
plt.colorbar(label=r"$\phi$")
plt.title("LBM D2Q9 - Solução ADR 2D")
plt.xlabel("X")
plt.ylabel("Y")
plt.tight_layout()
plt.savefig(f"{output}/solution.png", dpi=300, bbox_inches="tight")

# Plot Erro
plt.figure(figsize=(7, 6))
plt.imshow(erro_abs.T, origin="lower", extent=[0, Lx, 0, Ly], cmap="inferno")
plt.colorbar(label="Erro Absoluto")
plt.title("Erro Absoluto")
plt.tight_layout()
plt.savefig(f"{output}/error.png", dpi=300, bbox_inches="tight")

print(f"Resultados salvos em: {output}/")
