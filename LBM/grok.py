import numpy as np
import matplotlib.pyplot as plt
import os

# Pasta de saída
output_dir = "lbm_bubble_rising_frames"
os.makedirs(output_dir, exist_ok=True)

# Parâmetros
Nx, Ny = 100, 200           # Nx = largura (x), Ny = altura (y)
Nt = 1000
tau = 0.8
g = -0.00005                # ajuste se subir rápido demais

# Velocidades D2Q9
c = np.array([
    [0, 0],
    [1, 0], [0, 1], [-1, 0], [0, -1],
    [1, 1], [-1, 1], [-1, -1], [1, -1]
])
w = np.array([4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36])

# Distribuições e densidade inicial
f = np.zeros((9, Nx, Ny))
rho = np.ones((Nx, Ny))

# Bolha leve no fundo
radius = 20
center_x = Nx // 2
center_y = 40               # um pouco mais alto para dar espaço

# Correção aqui ↓
X, Y = np.meshgrid(np.arange(Nx), np.arange(Ny), indexing='ij')
mask = (X - center_x)**2 + (Y - center_y)**2 < radius**2
rho[mask] = 0.1

# Inicialização de equilíbrio (zero velocidade inicial)
for i in range(9):
    f[i, :, :] = w[i] * rho

# =============================================
# A partir daqui continua o loop temporal (exemplo mínimo)
# =============================================
for t in range(Nt):
    # Macroscópicas
    rho = np.sum(f, axis=0)
    ux = np.sum(c[:,0,None,None] * f, axis=0) / rho
    uy = np.sum(c[:,1,None,None] * f, axis=0) / rho

    # Gravidade (força de volume simples)
    uy += g

    # Equilíbrio
    feq = np.zeros_like(f)
    u2 = ux**2 + uy**2
    for i in range(9):
        cu = 3 * (c[i,0]*ux + c[i,1]*uy)
        feq[i] = w[i] * rho * (1 + cu + 0.5*cu**2 - 1.5*u2)

    # Colisão
    f -= (1/tau) * (f - feq)

    # Streaming (periódico)
    for i in range(9):
        f[i] = np.roll(f[i], shift=c[i,0], axis=0)
        f[i] = np.roll(f[i], shift=c[i,1], axis=1)

    # Salvar frames a cada 10 passos (menos arquivos)
    if t % 10 == 0:
        plt.figure(figsize=(5, 10))
        plt.imshow(rho.T, origin='lower', cmap='viridis')
        plt.colorbar(label='densidade')
        plt.title(f"t = {t:04d}  |  ρ_bolha ≈ {rho.min():.3f}")
        plt.tight_layout()
        plt.savefig(f"{output_dir}/frame_{t:04d}.png", dpi=100)
        plt.close()

print("Simulação terminada.")
