import numpy as np
import matplotlib.pyplot as plt
import os

# ======================
# Parâmetros
# ======================
nx, ny = 300, 100      # reduzi pra não travar seu PC
nt = 5000

tau = 0.6
ga = 1.111111111e-5

ux_in = 0.15
uy_in = 0.0
rho_in = 1.0

# D2Q9
ex = np.array([0,1,0,-1,0,1,-1,-1,1])
ey = np.array([0,0,1,0,-1,1,1,-1,-1])

w = np.array([4/9] + [1/9]*4 + [1/36]*4)

# ======================
# Campos
# ======================
rho = np.ones((ny, nx))
ux = np.zeros((ny, nx))
uy = np.zeros((ny, nx))

px = np.zeros((ny, nx))
py = np.zeros((ny, nx))

f = np.zeros((9, ny, nx))
feq = np.zeros_like(f)

nos_solidos = np.zeros((ny, nx))
nos_solidos[0,:] = 1
nos_solidos[-1,:] = 1

# ======================
# Funções
# ======================
def compute_feq(rho, ux, uy):
    feq = np.zeros((9, ny, nx))
    u2 = ux**2 + uy**2
    for k in range(9):
        cu = 3*(ex[k]*ux + ey[k]*uy)
        feq[k] = w[k]*rho*(1 + cu + 0.5*cu**2 - 1.5*u2)
    return feq


def macroscopic(f):
    rho = np.sum(f, axis=0)
    px = np.sum(f * ex[:,None,None], axis=0)
    py = np.sum(f * ey[:,None,None], axis=0)

    rho[rho < 1e-12] = 1e-12

    ux = px / rho
    uy = py / rho

    ux[nos_solidos == 1] = 0
    uy[nos_solidos == 1] = 0

    return rho, ux, uy


def entrada():
    global rho, ux, uy, f

    for i in range(1, ny-1):
        if nos_solidos[i,0] == 0:
            rho[i,0] = rho_in
            ux[i,0] = ux_in
            uy[i,0] = uy_in

            px = ux[i,0]*rho[i,0]
            py = uy[i,0]*rho[i,0]

            f[1,i,0] = f[3,i,0] + 2*px/3
            f[5,i,0] = f[7,i,0] + px/6 + py/2
            f[8,i,0] = f[6,i,0] + px/6 - py/2


def collision(f, feq):
    return f - (f - feq)/tau


def streaming(f):
    f_stream = np.zeros_like(f)

    for k in range(9):
        f_stream[k] = np.roll(f[k], shift=(ey[k], ex[k]), axis=(0,1))

    return f_stream


# ======================
# Inicialização
# ======================
feq = compute_feq(rho, ux, uy)
f = feq.copy()

# pasta de saída
os.makedirs("frames", exist_ok=True)

# ======================
# Loop principal
# ======================
for t in range(nt):

    # macroscópico
    rho, ux, uy = macroscopic(f)

    # entrada
    entrada()

    # força
    ux += tau * ga

    # equilíbrio
    feq = compute_feq(rho, ux, uy)

    # colisão
    f = collision(f, feq)

    # propagação
    f = streaming(f)

    # saída (direita)
    f[:, :, -1] = f[:, :, -2]

    # ======================
    # PLOT
    # ======================
    if t % 100 == 0:
        vel = np.sqrt(ux**2 + uy**2)

        plt.figure(figsize=(6,3))
        plt.imshow(vel, origin='lower')
        plt.colorbar(label='|u|')
        plt.title(f"Step {t}")
        plt.tight_layout()
        plt.savefig(f"frames/frame_{t:05d}.png")
        plt.close()

        print(f"Step {t}, max u = {vel.max():.5f}")
