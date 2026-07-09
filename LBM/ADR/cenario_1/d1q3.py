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
with open("d1q3.yaml", "r") as file:
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
# DOMÍNIO 1D
# ==========================================================
Nx = config["dominio"]["Nx"]
Lx = config["dominio"]["Lx"]
dx = Lx / (Nx - 1)
x = np.linspace(0, Lx, Nx)

# ==========================================================
# PARÂMETROS FÍSICOS
# ==========================================================
ux = config["fisica"]["ux"]
tau = config["fisica"]["tau"]
k_reaction = config["fisica"].get("k_reaction", 0.0)
cs2 = 1.0 / 3.0
D = cs2 * (tau - 0.5)

# ==========================================================
# MANUFACTURED SOLUTION (MMS)
# ==========================================================
kx = 2 * np.pi / Lx
phi_exact = np.sin(kx * x)

dphi_dx = kx * np.cos(kx * x)
d2phi_dx2 = - (kx**2) * phi_exact

S = ux * dphi_dx - D * d2phi_dx2 + k_reaction * phi_exact

ux_lattice = ux * dx
S_lattice = S * dx**2

print("Configuração carregada do d1q3.yaml:")
print(f"Nx = {Nx} | Lx = {Lx}")
print(f"ux = {ux} | tau = {tau} | D = {D:.6f} | k_reaction = {k_reaction}")
print(f"Fonte máxima = {np.max(np.abs(S)):.4f}")

# ==========================================================
# LATTICE D1Q3
# ==========================================================
c = np.array([-1, 0, 1])
w = np.array([1/6, 2/3, 1/6])
Q = 3

# Inicialização
phi = phi_exact.copy() * 0.01
f = np.zeros((Q, Nx))

def get_equilibrium(phi_field):
    feq = np.zeros((Q, Nx))
    for i in range(Q):
        cu = c[i] * ux_lattice / cs2
        feq[i] = w[i] * phi_field * (1.0 + cu)
    return feq

f = get_equilibrium(phi)

# ==========================================================
# SIMULAÇÃO
# ==========================================================
tol = float(config["simulacao"]["tol"])        # ← CORREÇÃO AQUI
maxIter = int(config["simulacao"]["maxIter"])
history = []

start = time.perf_counter()

for it in tqdm(range(maxIter)):
    phi_old = np.sum(f, axis=0) + 0.5 * S_lattice
    
    feq = get_equilibrium(phi_old)
    for i in range(Q):
        f[i] = f[i] - (f[i] - feq[i]) / tau + w[i] * (1.0 - 0.5 / tau) * S_lattice
    
    # Streaming
    for i in range(Q):
        f[i] = np.roll(f[i], c[i])
    
    # Dirichlet BC
    for i in range(Q):
        cu = c[i] * ux_lattice / cs2
        f[i, 0]  = w[i] * phi_exact[0]  * (1.0 + cu)
        f[i, -1] = w[i] * phi_exact[-1] * (1.0 + cu)
    
    phi = np.sum(f, axis=0) + 0.5 * S_lattice
    
    err = np.sqrt(np.sum((phi - phi_old)**2) / (np.sum(phi_old**2) + 1e-15))
    history.append(err)
    
    if err < tol:
        break

end = time.perf_counter()

# ==========================================================
# PÓS-PROCESSAMENTO
# ==========================================================
erro_abs = np.abs(phi - phi_exact)
L2_rel = np.sqrt(np.sum((phi - phi_exact)**2) / np.sum(phi_exact**2))

print("\n" + "="*70)
print("LBM D1Q3 - ADR 1D Steady State (MMS)")
print("="*70)
print(f"Iterações : {it}")
print(f"Tempo     : {end - start:.3f} s")
print(f"L2 relativo: {L2_rel:.2e}")
print("="*70)

# Salvar resultados
with open(f"{output}/solution_1d.csv", "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["x", "phi_LBM", "phi_exact", "error"])
    for i in range(Nx):
        writer.writerow([x[i], float(phi[i]), float(phi_exact[i]), float(erro_abs[i])])

# Gráficos
plt.figure(figsize=(10, 6))
plt.plot(x, phi_exact, 'b-', linewidth=2.5, label="Solução Exata")
plt.plot(x, phi, 'r--', linewidth=2, label="LBM D1Q3")
plt.title("LBM D1Q3 vs Solução Exata\nEquação ADR 1D Steady-State")
plt.xlabel("Posição x")
plt.ylabel(r"$\phi(x)$")
plt.grid(True, alpha=0.5)
plt.legend()
save_fig("comparacao_lbm_exact.png")

plt.figure(figsize=(10, 4))
plt.plot(x, erro_abs, 'g-')
plt.title("Erro Absoluto")
plt.xlabel("x")
plt.ylabel("Erro")
plt.grid(True)
save_fig("erro.png")

plt.figure(figsize=(8, 4))
plt.semilogy(history)
plt.title("Convergência")
plt.xlabel("Iterações")
plt.ylabel("Resíduo")
plt.grid(True)
save_fig("convergencia.png")

print(f"Resultados salvos em: {output}/")
