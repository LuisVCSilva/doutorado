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
    plt.savefig(
        f"{output}/{filename}",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()


# ==========================================================
# DOMÍNIO
# ==========================================================

Nx = config["dominio"]["Nx"]
Ny = config["dominio"]["Ny"]

Lx = config["dominio"]["Lx"]
Ly = config["dominio"]["Ly"]

dx = Lx/(Nx-1)
dy = Ly/(Ny-1)

x = np.linspace(0,Lx,Nx)
y = np.linspace(0,Ly,Ny)

X,Y = np.meshgrid(
    x,
    y,
    indexing="ij"
)


# ==========================================================
# FÍSICA
# ==========================================================

ux = config["fisica"]["ux"]
uy = config["fisica"]["uy"]

tau = config["fisica"]["tau"]

cs2 = 1/3


D = cs2*(tau-0.5)


Pe = (
    np.sqrt(ux**2 + uy**2)
    * Lx
    / D
)


# ==========================================================
# LATTICE
# ==========================================================

ux_lattice = ux
uy_lattice = uy

D_lattice = D


# ==========================================================
# MMS
# ==========================================================

X_lat,Y_lat = np.meshgrid(
    np.arange(Nx),
    np.arange(Ny),
    indexing="ij"
)


kx = np.pi/(Nx-1)
ky = np.pi/(Ny-1)


phi_exact = (
    np.sin(kx*X_lat)
    *
    np.sin(ky*Y_lat)
)


dphi_dx = (
    kx
    *
    np.cos(kx*X_lat)
    *
    np.sin(ky*Y_lat)
)


dphi_dy = (
    ky
    *
    np.sin(kx*X_lat)
    *
    np.cos(ky*Y_lat)
)


laplacian = (
    -(kx**2+ky**2)
    *
    phi_exact
)


# Fonte MMS

S_lattice = (
    ux_lattice*dphi_dx
    +
    uy_lattice*dphi_dy
    -
    D_lattice*laplacian
)



# ==========================================================
# SIMULAÇÃO
# ==========================================================

tol = config["simulacao"]["tol"]
maxIter = config["simulacao"]["maxIter"]


print("-----------------------")
print("Configuração")
print("-----------------------")
print(f"Nx Ny = {Nx} {Ny}")
print(f"D     = {D}")
print(f"Pe    = {Pe}")
print(f"tau   = {tau}")
print(f"maxIt = {maxIter}")


# ==========================================================
# D2Q5
# ==========================================================

c = np.array([
    [0,0],
    [1,0],
    [-1,0],
    [0,1],
    [0,-1]
])


w = np.array([
    1/3,
    1/6,
    1/6,
    1/6,
    1/6
])


Q = 5


# ==========================================================
# EQUILÍBRIO D2Q5 PARA ADR
# ==========================================================

def equilibrium(phi):

    feq = np.zeros(
        (Q, Nx, Ny)
    )

    for i in range(Q):

        cu = (
            c[i,0]*ux_lattice
            +
            c[i,1]*uy_lattice
        )

        feq[i] = (
            w[i]
            *
            phi
            *
            (
                1.0
                +
                cu/cs2
            )
        )

    return feq



# ==========================================================
# INICIALIZAÇÃO
# ==========================================================

phi = phi_exact.copy()

f = equilibrium(phi)

history = []

start = time.perf_counter()



# ==========================================================
# LOOP LBM
# ==========================================================

for it in tqdm(range(maxIter)):


    # variável macroscópica

    phi_old = (
        np.sum(f,axis=0)
        +
        0.5*S_lattice
    )


    # equilíbrio

    feq = equilibrium(phi_old)



    # ------------------------------------------------------
    # Colisão BGK + Fonte de Guo
    # ------------------------------------------------------

    for i in range(Q):

        f[i] = (
            f[i]
            -
            (f[i]-feq[i])/tau
            +
            w[i]
            *
            (
                1.0-0.5/tau
            )
            *
            S_lattice
        )



    # ------------------------------------------------------
    # Streaming
    # ------------------------------------------------------

    for i in range(Q):

        f[i] = np.roll(
            f[i],
            c[i,0],
            axis=0
        )

        f[i] = np.roll(
            f[i],
            c[i,1],
            axis=1
        )



    # ------------------------------------------------------
    # Condição de contorno Dirichlet
    # usando equilíbrio completo
    # ------------------------------------------------------

    feq_bc = equilibrium(phi_exact)


    for i in range(Q):

        # esquerda
        f[i,0,:] = feq_bc[i,0,:]

        # direita
        f[i,-1,:] = feq_bc[i,-1,:]

        # inferior
        f[i,:,0] = feq_bc[i,:,0]

        # superior
        f[i,:,-1] = feq_bc[i,:,-1]



    # ------------------------------------------------------
    # Atualização macroscópica
    # ------------------------------------------------------

    phi = (
        np.sum(f,axis=0)
        +
        0.5*S_lattice
    )



    # ------------------------------------------------------
    # Resíduo
    # ------------------------------------------------------

    err = np.sqrt(
        np.sum((phi-phi_old)**2)
        /
        (
            np.sum(phi_old**2)
            +
            1e-15
        )
    )


    history.append(err)


    if err < tol:
        break



end = time.perf_counter()



# ==========================================================
# ERRO MMS
# ==========================================================

erro = np.abs(
    phi-phi_exact
)


L2 = np.sqrt(
    np.sum((phi-phi_exact)**2)
    /
    np.sum(phi_exact**2)
)



print("\n")
print("==============================")
print("LBM D2Q5 ADR 2D")
print("==============================")
print(f"Iterações : {it}")
print(f"Tempo      : {end-start:.4f}")
print(f"D          : {D:.6e}")
print(f"Pe         : {Pe:.6f}")
print(f"L2         : {L2:.6e}")
print("==============================")


# ==========================================================
# CSV
# ==========================================================

with open(
    output+"/metrics.csv",
    "w",
    newline=""
) as file:

    writer=csv.writer(file)

    writer.writerow([
        "Nx",
        "Ny",
        "tau",
        "D",
        "Pe",
        "iterations",
        "time",
        "L2"
    ])

    writer.writerow([
        Nx,
        Ny,
        tau,
        D,
        Pe,
        it,
        end-start,
        L2
    ])



# ==========================================================
# FIGURAS
# ==========================================================


plt.figure(figsize=(7,6))

plt.imshow(
    phi.T,
    origin="lower",
    extent=[0,Lx,0,Ly]
)

plt.colorbar()

plt.title(
    "LBM D2Q5"
)

save_fig(
    "solution.png"
)



plt.figure(figsize=(7,6))

plt.imshow(
    erro.T,
    origin="lower",
    extent=[0,Lx,0,Ly]
)

plt.colorbar()

plt.title(
    "Erro absoluto"
)

save_fig(
    "error.png"
)



plt.figure(figsize=(7,4))

plt.semilogy(history)

plt.xlabel("Iteração")
plt.ylabel("Resíduo")

plt.grid()

save_fig(
    "convergence.png"
)



print(
    f"Resultados em {output}"
)
