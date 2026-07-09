import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import matplotlib
matplotlib.use("Agg")

def save_fig(filename):
    plt.tight_layout()
    plt.savefig(f"{output}/{filename}", dpi=300, bbox_inches="tight")
    plt.close()

# ==========================================================
# Configuração do diretório de saída
output = "dados/solucao_exata_2d"
os.makedirs(output, exist_ok=True)

# ==========================================================
# Configuração da Malha Espacial (Continuum)
Nx = 100
Ny = 100
Lx = 1.0
Ly = 1.0
x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)
X, Y = np.meshgrid(x, y, indexing="ij")

# ==========================================================
# Cálculo Analítico da Solução e Suas Derivadas
# Solução Exata Alvo: phi(x,y) = sin(2*pi*x) * sin(2*pi*y)
phi_exact = np.sin(2 * np.pi * X / Lx) * np.sin(2 * np.pi * Y / Ly)

# Primeira derivada em relação a X
dphidx = (2 * np.pi / Lx) * np.cos(2 * np.pi * X / Lx) * np.sin(2 * np.pi * Y / Ly)

# Primeira derivada em relação a Y
dphidy = (2 * np.pi / Ly) * np.sin(2 * np.pi * X / Lx) * np.cos(2 * np.pi * Y / Ly)

# Laplaciano: d2phi/dx2 + d2phi/dy2
laplacian = -((2 * np.pi / Lx)**2 + (2 * np.pi / Ly)**2) * phi_exact

# ==========================================================
# Salvar os Dados Numéricos em um Arquivo CSV
with open(output + "/exact_fields.csv", "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["x", "y", "phi_exact", "dphi_dx", "dphi_dy", "laplacian"])
    for i in range(Nx):
        for j in range(Ny):
            writer.writerow([X[i,j], Y[i,j], phi_exact[i,j], dphidx[i,j], dphidy[i,j], laplacian[i,j]])

# ==========================================================
# Geração dos Gráficos dos Campos Analíticos

# 1. Plot da Solução Exata
plt.figure(figsize=(7,6))
plt.imshow(phi_exact.T, origin="lower", extent=[0, Lx, 0, Ly], cmap="viridis")
plt.colorbar(label=r"$\phi_{exata}$")
plt.title(r"Solução Exata Analítica $\phi(x,y)$")
plt.xlabel("X")
plt.ylabel("Y")
save_fig("phi_exact.png")

# 2. Plot do Gradiente em X
plt.figure(figsize=(7,6))
plt.imshow(dphidx.T, origin="lower", extent=[0, Lx, 0, Ly], cmap="RdBu_r")
plt.colorbar(label=r"$\partial\phi/\partial x$")
plt.title(r"Gradiente Espacial $\frac{\partial\phi}{\partial x}$")
plt.xlabel("X")
plt.ylabel("Y")
save_fig("gradient_x.png")

# 3. Plot do Laplaciano
plt.figure(figsize=(7,6))
plt.imshow(laplacian.T, origin="lower", extent=[0, Lx, 0, Ly], cmap="coolwarm")
plt.colorbar(label=r"$\nabla^2\phi$")
plt.title(r"Laplaciano Analítico $\nabla^2\phi$")
plt.xlabel("X")
plt.ylabel("Y")
save_fig("laplacian.png")

print("="*50)
print("Gerador de Solução Exata 2D")
print("="*50)
print(f"Malha gerada       : {Nx} x {Ny}")
print(f"Resultados salvos em: {output}/")
print("="*50)
