import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import yaml
import matplotlib
matplotlib.use("Agg")

def load_config(config_file):
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def save_fig(filename, output):
    plt.tight_layout()
    plt.savefig(f"{output}/{filename}", dpi=300, bbox_inches="tight")
    plt.close()

# ==========================================================
# CARREGAR CONFIGURAÇÃO
# ==========================================================
config = load_config("fdm1d.yaml")
Nx = config['dominio']['Nx']
L = config['dominio']['L']
u = config['parametros']['u']
D = config['parametros']['D']
R = config['parametros']['R']
output = config['saida']['diretorio']

os.makedirs(output, exist_ok=True)

dx = L / (Nx - 1)
x = np.linspace(0, L, Nx)

# ==========================================================
# SOLUÇÃO ANALÍTICA (MMS)
# ==========================================================
phi_exact = np.sin(2 * np.pi * x / L)
dphi_exact = (2 * np.pi / L) * np.cos(2 * np.pi * x / L)
d2phi_exact = -(2 * np.pi / L)**2 * np.sin(2 * np.pi * x / L)

# Fonte para Manufactured Solution
S = u * dphi_exact - D * d2phi_exact + R * phi_exact

# ==========================================================
# MONTAGEM DO SISTEMA FDM (2ª ordem central)
# ==========================================================
A = np.zeros((Nx, Nx))
b = np.zeros(Nx)

# Condições de contorno Dirichlet
A[0, 0] = 1.0
b[0] = phi_exact[0]
A[-1, -1] = 1.0
b[-1] = phi_exact[-1]

for i in range(1, Nx - 1):
    A[i, i-1] = D / dx**2 + u / (2 * dx)
    A[i, i]   = -2 * D / dx**2 + R
    A[i, i+1] = D / dx**2 - u / (2 * dx)
    b[i] = -S[i]

# Solução do sistema linear
phi = np.linalg.solve(A, b)

# ==========================================================
# PÓS-PROCESSAMENTO
# ==========================================================
erro_abs = np.abs(phi - phi_exact)
L2_rel = np.sqrt(np.sum((phi - phi_exact)**2) / np.sum(phi_exact**2))

print("="*60)
print("FDM 1D - ADR Steady State")
print("="*60)
print(f"L2 relativo: {L2_rel:.2e}")
print(f"Resultados salvos em: {output}/")
print("="*60)

# Salvar CSV
with open(f"{output}/solution.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["x", "phi_fdm", "phi_exact", "error_abs"])
    for i in range(Nx):
        writer.writerow([x[i], phi[i], phi_exact[i], erro_abs[i]])

# ==========================================================
# GRÁFICOS (mesmo estilo do LBM)
# ==========================================================

# 1. Comparação LBM-style
plt.figure(figsize=(10, 6))
plt.plot(x, phi_exact, 'b-', linewidth=2.5, label="Solução Exata")
plt.plot(x, phi, 'r--', linewidth=2, label="FDM")
plt.title("FDM vs Solução Exata\nEquação ADR 1D Steady-State")
plt.xlabel("Posição x")
plt.ylabel(r"$\phi(x)$")
plt.grid(True, alpha=0.5)
plt.legend()
save_fig("comparacao_fdm_exact.png", output)

# 2. Erro Absoluto
plt.figure(figsize=(10, 4))
plt.plot(x, erro_abs, 'g-')
plt.title("Erro Absoluto - FDM")
plt.xlabel("x")
plt.ylabel("Erro")
plt.grid(True)
save_fig("erro.png", output)

# 3. (Opcional) Solução simples
plt.figure(figsize=(8, 5))
plt.plot(x, phi, label="FDM", color='red')
plt.plot(x, phi_exact, "--", label="Analítica", color='blue')
plt.xlabel("x")
plt.ylabel(r"$\phi$")
plt.legend()
plt.grid(True)
save_fig("solution.png", output)

print("Todos os plots gerados com sucesso!")
