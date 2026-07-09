import numpy as np
import matplotlib.pyplot as plt
import csv
import time
import os
import yaml
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import matplotlib
matplotlib.use("Agg")

def load_config(config_file):
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def save_fig(filename, output):
    plt.tight_layout()
    plt.savefig(f"{output}/{filename}", dpi=300, bbox_inches="tight")
    plt.close()

# Carregar configurações
config = load_config("fdm2d.yaml")
Nx, Ny = config['dominio']['Nx'], config['dominio']['Ny']
Lx, Ly = config['dominio']['Lx'], config['dominio']['Ly']
ux, uy = config['fisica']['ux'], config['fisica']['uy']
print(config['fisica'])
D = float(config['fisica']['D'])
R = float(config['fisica']['R'])
output = config['saida']['diretorio']

os.makedirs(output, exist_ok=True)
dx, dy = Lx/(Nx-1), Ly/(Ny-1)
x, y = np.linspace(0, Lx, Nx), np.linspace(0, Ly, Ny)
X, Y = np.meshgrid(x, y, indexing="ij")
N = Nx * Ny



# ==========================================================
# PARÂMETROS ADR
# ==========================================================

ux = 0.10
uy = 0.10


D = 0.01


R = -1.0


Pe = np.sqrt(
    ux**2+uy**2
)*Lx/D



# ==========================================================
# SOLUÇÃO ANALÍTICA
# ==========================================================


phi_exact = (

    np.sin(2*np.pi*X/Lx)
    *
    np.sin(2*np.pi*Y/Ly)

)



# derivadas

dphidx = (

    2*np.pi/Lx
    *
    np.cos(2*np.pi*X/Lx)
    *
    np.sin(2*np.pi*Y/Ly)

)



dphidy = (

    2*np.pi/Ly
    *
    np.sin(2*np.pi*X/Lx)
    *
    np.cos(2*np.pi*Y/Ly)

)



lap = (

    -(2*np.pi/Lx)**2*phi_exact

    -

    (2*np.pi/Ly)**2*phi_exact

)



# Fonte manufaturada

S = (
    ux * dphidx +
    uy * dphidy -
    D * lap -
    R * phi_exact
)

# Perturbação: Dois pulsos opostos fortes
pulse = np.zeros((Nx, Ny))
pulse[Nx//3, Ny//3] = 1000.0      # pulso positivo
pulse[2*Nx//3, 2*Ny//3] = -1000.0 # pulso negativo

S = S + pulse

# ==========================================================
# MATRIZ FDM
# ==========================================================


A = lil_matrix(
    (N,N)
)


b = np.zeros(N)



def idx(i,j):

    return i*Ny+j




# ==========================================================
# MONTAGEM
# ==========================================================


for i in range(Nx):

    for j in range(Ny):


        k = idx(i,j)



        # --------------------------
        # fronteira Dirichlet
        # --------------------------

        if (

            i==0
            or i==Nx-1
            or j==0
            or j==Ny-1

        ):


            A[k,k]=1.0

            b[k]=phi_exact[i,j]



        else:


            # coeficientes difusão

            aw = D/dx**2 + ux/(2*dx)

            ae = D/dx**2 - ux/(2*dx)


            an = D/dy**2 + uy/(2*dy)

            ass = D/dy**2 - uy/(2*dy)



            ap = (

                -2*D/dx**2

                -

                2*D/dy**2

                +

                R

            )



            A[k,idx(i-1,j)] = aw

            A[k,idx(i+1,j)] = ae


            A[k,idx(i,j-1)] = ass

            A[k,idx(i,j+1)] = an


            A[k,k]=ap



            b[k]=-S[i,j]




# ==========================================================
# SOLUÇÃO
# ==========================================================


start=time.perf_counter()


phi_vec = spsolve(
    A.tocsr(),
    b
)


end=time.perf_counter()



phi = phi_vec.reshape(
    Nx,
    Ny
)



# ==========================================================
# ERRO
# ==========================================================


erro=np.abs(
    phi-phi_exact
)



L2=np.sqrt(

    np.sum(
        (phi-phi_exact)**2
    )

    /

    np.sum(phi_exact**2)

)



print()

print("==============================")
print("FDM ADR 2D")
print("==============================")
print("Nx:",Nx)
print("Ny:",Ny)
print("Tempo:",end-start)
print("Pe:",Pe)
print("L2:",L2)
print("==============================")



# ==========================================================
# CSV MÉTRICAS
# ==========================================================

with open(
    output+"/metrics.csv",
    "w",
    newline=""
) as f:


    writer=csv.writer(f)


    writer.writerow([

        "Nx",
        "Ny",
        "D",
        "ux",
        "uy",
        "R",
        "Pe",
        "time",
        "L2"

    ])


    writer.writerow([

        Nx,
        Ny,
        D,
        ux,
        uy,
        R,
        Pe,
        end-start,
        L2

    ])




# ==========================================================
# CSV SOLUÇÃO
# ==========================================================


with open(
    output+"/solution.csv",
    "w",
    newline=""
) as f:


    writer=csv.writer(f)


    writer.writerow([

        "x",
        "y",
        "phi_FDM",
        "phi_exact",
        "error"

    ])



    for i in range(Nx):

        for j in range(Ny):

            writer.writerow([

                X[i,j],
                Y[i,j],
                phi[i,j],
                phi_exact[i,j],
                erro[i,j]

            ])




# ==========================================================
# SOLUÇÃO FDM
# ==========================================================


plt.figure(figsize=(7,6))

plt.imshow(
    phi.T,
    origin="lower",
    extent=[
        0,Lx,
        0,Ly
    ],
    aspect="auto"
)


plt.colorbar()

plt.title(
    "ADR 2D Steady State - FDM"
)


save_fig(
    "solution.png",output
)




# ==========================================================
# SOLUÇÃO ANALÍTICA
# ==========================================================


plt.figure(figsize=(7,6))


plt.imshow(
    phi_exact.T,
    origin="lower",
    extent=[
        0,Lx,
        0,Ly
    ],
    aspect="auto"
)


plt.colorbar()


plt.title(
    "Solução Analítica"
)


save_fig(
    "analytic.png",output
)




# ==========================================================
# ERRO
# ==========================================================


plt.figure(figsize=(7,6))


plt.imshow(
    erro.T,
    origin="lower",
    extent=[
        0,Lx,
        0,Ly
    ],
    aspect="auto"
)


plt.colorbar()


plt.title(
    "Erro absoluto"
)


save_fig(
    "error.png",output
)

