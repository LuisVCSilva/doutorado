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
output = "dados/solucao_exata_1d"
os.makedirs(output, exist_ok=True)

# ==========================================================
# Configuração da Malha Espacial 1D (Continuum)
Nx = 100
Lx = 1.0
x = np.linspace(0, Lx, Nx)

# ==========================================================
# Cálculo Analítico da Solução e Suas Derivadas
# Solução Exata Alvo: phi(x) = sin(2*pi*x)
phi_exact = np.sin(2 * np.pi * x / Lx)

# Primeira derivada (Gradiente 1D)
dphidx = (2 * np.pi / Lx) * np.cos(2 * np.pi * x / Lx)

# Segunda derivada (Laplaciano 1D)
d2phidx2 = -((2 * np.pi / Lx)**2) * phi_exact

# ==========================================================
# Salvar os Dados Numéricos em um Arquivo CSV
with open(output + "/exact_fields_1d.csv", "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["x", "phi_exact", "dphi_dx", "d2phi_dx2"])
    for i in range(Nx):
        writer.writerow([x[i], phi_exact[i], dphidx[i], d2phidx2[i]])

# ==========================================================
# Geração dos Gráficos Analíticos (Curvas 1D)

# Gráfico Combinado para visualização clara das relações moleculares
plt.figure(figsize=(8, 5))
plt.plot(x, phi_exact, label=r"Solução Exata $\phi(x)$", color="darkblue", linewidth=2)
plt.plot(x, dphidx, label=r"1ª Derivada $\frac{d\phi}{dx}$", color="crimson", linestyle="--")
plt.plot(x, d2phidx2, label=r"2ª Derivada $\frac{d^2\phi}{dx^2}$", color="forestgreen", linestyle=":")

plt.title("Campos Analíticos 1D Exatos")
plt.xlabel("Posição (x)")
plt.ylabel("Amplitude")
plt.grid(True, which="both", linestyle="--", alpha=0.5)
plt.legend(loc="upper right")
save_fig("exact_solution_1d.png")

print("="*50)
print("Gerador de Solução Exata 1D")
print("="*50)
print(f"Malha gerada       : {Nx} nós")
print(f"Resultados salvos em: {output}/")
print("="*50)
