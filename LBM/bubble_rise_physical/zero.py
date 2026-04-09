import numpy as np
import matplotlib.pyplot as plt
import os

# =============================
# pasta
# =============================
output_dir = "bubble_clean"
os.makedirs(output_dir, exist_ok=True)

# =============================
# domínio
# =============================
Nx, Ny = 180, 300
Nt = 3000

tau = 0.65

# =============================
# D2Q9
# =============================
c = np.array([
    [0,0],
    [1,0],[0,1],[-1,0],[0,-1],
    [1,1],[-1,1],[-1,-1],[1,-1]
])

w = np.array([
    4/9,
    1/9,1/9,1/9,1/9,
    1/36,1/36,1/36,1/36
])

opp = np.array([0,3,4,1,2,7,8,5,6])

# =============================
# fluido
# =============================
rho = np.ones((Nx,Ny))
ux = np.zeros((Nx,Ny))
uy = np.zeros((Nx,Ny))

# =============================
# bolha (campo fase)
# =============================
phi = np.ones((Nx,Ny))

cx, cy = Nx//2, Ny//6
R = 18

for x in range(Nx):
    for y in range(Ny):
        r = np.sqrt((x-cx)**2 + (y-cy)**2)
        if r < R:
            phi[x,y] = 0.0

# =============================
# distribuição inicial
# =============================
f = np.zeros((9,Nx,Ny))

for i in range(9):
    f[i] = w[i]*rho

# =============================
# loop temporal
# =============================
for t in range(Nt):

    # -------------------------
    # macroscópicas
    # -------------------------
    rho = np.sum(f, axis=0)

    ux = np.sum(c[:,0,None,None]*f, axis=0)/rho
    uy = np.sum(c[:,1,None,None]*f, axis=0)/rho

    # -------------------------
    # empuxo local da bolha
    # -------------------------
    buoyancy = 2e-5*(1-phi)

    uy += buoyancy

    # -------------------------
    # deformação da bolha
    # -------------------------
    gradx = np.roll(phi,-1,axis=0)-np.roll(phi,1,axis=0)
    grady = np.roll(phi,-1,axis=1)-np.roll(phi,1,axis=1)

    ux += 0.002*gradx
    uy += 0.002*grady

    # -------------------------
    # equilíbrio
    # -------------------------
    feq = np.zeros_like(f)

    u2 = ux**2 + uy**2

    for i in range(9):
        cu = 3*(c[i,0]*ux + c[i,1]*uy)
        feq[i] = w[i]*rho*(1 + cu + 0.5*cu**2 - 1.5*u2)

    # -------------------------
    # colisão
    # -------------------------
    f = f - (1/tau)*(f-feq)

    # -------------------------
    # streaming
    # -------------------------
    for i in range(9):
        f[i] = np.roll(np.roll(f[i], c[i,0], axis=0), c[i,1], axis=1)

    # -------------------------
    # bounce-back fundo
    # -------------------------
    for i in range(9):
        f[i,:,0] = f[opp[i],:,0]

    # -------------------------
    # topo aberto
    # -------------------------
    for i in range(9):
        f[i,:,-1] = f[i,:,-2]

    # -------------------------
    # advecção da bolha
    # -------------------------
    phi_new = phi.copy()

    for x in range(1,Nx-1):
        for y in range(1,Ny-1):
            dx = -ux[x,y]*(phi[x+1,y]-phi[x-1,y])*0.5
            dy = -uy[x,y]*(phi[x,y+1]-phi[x,y-1])*0.5

            phi_new[x,y] += 0.2*(dx+dy)

    phi = np.clip(phi_new,0,1)

    # -------------------------
    # suavização interface
    # -------------------------
    phi = 0.96*phi + 0.01*(
        np.roll(phi,1,0)+np.roll(phi,-1,0)+
        np.roll(phi,1,1)+np.roll(phi,-1,1)
    )

    # -------------------------
    # salvar
    # -------------------------
    if t % 20 == 0:
        plt.figure(figsize=(6,10))
        plt.imshow(phi.T, origin='lower')
        plt.colorbar(label='phase')
        plt.title(f"Bubble rise t={t}")
        plt.savefig(f"{output_dir}/frame_{t:04d}.png")
        plt.close()
