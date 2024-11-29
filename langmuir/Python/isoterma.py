import sys, json, time, os,subprocess, random, math
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
global k
k = 0

def solucao_analitica(t, rA, rD, rc):
    return (rA / (rA + rD)) * (1 - np.exp(-(rA + rD) * t)) + rc * np.exp(-(rA + rD) * t)
    
    
def inicializar_rede(Lx, Ly, rc):
    isite = (np.random.rand(Lx, Ly) < rc).astype(bool)
    ic = np.sum(isite)
    return isite, ic
    

def simular(Lx, Ly, LT, rA, rD, isite, ic, intervalo_snapshot, nome_arquivo_saida, modo=["csv", "plot"]):
    t = 0.0
    total_iteracoes = int(LT / (1 / (ic * rD + (Lx * Ly - ic) * rA))) 
    historico_theta = np.zeros(total_iteracoes)
    historico_tempo = np.linspace(0, LT, total_iteracoes)
    contagem = 0
    valores_aleatorios = np.random.rand(total_iteracoes) 
    cond_d = rD / (rA + rD)
    cond_a = rA / (rA + rD)
    area = (Lx * Ly)

    indices_i = np.random.randint(0, Lx, total_iteracoes)
    indices_j = np.random.randint(0, Ly, total_iteracoes)
    valores_aleatorios = np.random.rand(total_iteracoes) 

    for _ in tqdm(range(total_iteracoes), desc=""):
        i, j = indices_i[_], indices_j[_]
        r = valores_aleatorios[_]
        
        if isite[i, j] != 0: 
            if r <= (cond_d): 
                isite[i, j] = 0
                ic -= 1
                dt = -math.log(r) / (ic * rD + (Lx * Ly - ic) * rA)
                t += dt
        else:  
            if r <= (cond_a):
                isite[i, j] = 1
                ic += 1
                dt = -math.log(r) / (ic * rD + (Lx * Ly - ic) * rA)
                t += dt
                        
        rc = ic / area  
        historico_theta[contagem] = rc 
        historico_tempo[contagem] = t     
        contagem = contagem + 1

        if _ % intervalo_snapshot == 0:
            plotar_combinado(isite, t, historico_tempo[0:_], historico_theta[0:_], rc, nome_arquivo_saida, modo)

        if t >= LT:  # Parar se o tempo máximo for atingido
            historico_theta = historico_theta[:contagem]
            historico_tempo = historico_tempo[:contagem]
            plotar_combinado(isite, t, historico_tempo[0:_], historico_theta[0:_], rc, nome_arquivo_saida, modo)
            break

    return historico_tempo, historico_theta

def plotar_combinado(snapshot, t, historico_tempo, historico_theta, rc, nome_arquivo, modo=["csv", "plot"]):
    global k
    k += 1
    if "csv" in modo:
        historico_df = pd.DataFrame({'Tempo': historico_tempo, 'Theta': historico_theta})

        arquivo_historico = f"{nome_arquivo}/time_series/theta.csv"
        if os.path.exists(arquivo_historico):
            historico_df.to_csv(arquivo_historico, mode='a', header=False, index=False) 
        else:
            historico_df.to_csv(arquivo_historico, mode='w', header=True, index=False) 
        

        rede_df = pd.DataFrame(snapshot)  # Criar um DataFrame a partir do snapshot
        rede_df.to_csv(f'{nome_arquivo}/lattice/estado_rede_{k}.csv', index=False, header=False) 
    if "plot" in modo:

        fig, axs = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True) 
    

        axs[0].imshow(snapshot, cmap='binary', origin='lower', interpolation='nearest')
        axs[0].set_title(f"Autômatos de gás -- ${Lx}x{Ly}$ $t = {t:.2f} \\ s$", fontsize=14)
        axs[0].axis('off')

        theta_analitico = solucao_analitica(historico_tempo, rA, rD, rc)
        ultimo_theta_analitico = float('nan') if len(theta_analitico)==0 else theta_analitico[-1]

        axs[1].plot(historico_tempo, theta_analitico, label="Solução Analítica", linestyle='--', color="red")

        axs[1].plot(historico_tempo, historico_theta, label='Cobertura', color='blue')

        axs[1].set_title(f'KMC Langmuir  $T_a = {rA:.2f} \\ T_d = {rD:.2f} \\ MC(r_c) = {rc * 100:.2f}\\% \\ Exato(r_c) = {ultimo_theta_analitico * 100:.2f}$')
        axs[1].set_xlabel('Tempo t')
        axs[1].set_ylabel('Cobertura (theta)')
        axs[1].grid()
    

        axs[1].legend(loc='upper right') 
    

        arquivo_rede = f"{nome_arquivo}/plot/snapshot_combinado_{k}.png"
        plt.savefig(arquivo_rede, bbox_inches='tight', pad_inches=0.1)
        plt.close(fig)  
def plotar_resultados(historico_tempo, historico_theta, rA, rD, rc):
    plt.figure(figsize=(10, 5))


    plt.plot(historico_tempo, historico_theta, label="Dados Simulados", color="blue")


    array_tempo = np.linspace(0, historico_tempo[-1], len(historico_tempo))
    theta_analitico = solucao_analitica(array_tempo, rA, rD, rc)


    plt.plot(array_tempo, theta_analitico, label="Solução Analítica", linestyle='--', color="red")


    plt.title('Cobertura ao Longo do Tempo')
    plt.xlabel('Tempo')
    plt.ylabel('Cobertura (theta)')
    plt.legend()
    plt.grid()
    plt.show()

def main(nome_arquivo_config):

    with open(nome_arquivo_config, 'r') as f:
        config = json.load(f)


    global Lx, Ly, LT, rc, rA, rD
    Lx = config.get("Lx", 2**6)
    Ly = config.get("Ly", 2**6)
    LT = config.get("LT", int(1e+3))
    rc = config.get("rc", 0.1)
    rA = config.get("rA", 2.0)
    rD = config.get("rD", 1.0)
    modo = config.get("modo", "csv,plot")
    intervalo_snapshot = config.get("snapshot_intervalo", 1000)
    # Configuração do diretório para saída
    timestamp = int(time.time())
    dirname = "sim_" + str(timestamp)
    os.makedirs(dirname, exist_ok=True)
    os.makedirs(dirname+"/plot", exist_ok=True)
    os.makedirs(dirname+"/lattice", exist_ok=True)
    os.makedirs(dirname+"/time_series", exist_ok=True)
    os.makedirs(dirname+"/video", exist_ok=True)

    with open(f"{dirname}/params.json", 'w') as json_file:
        json.dump(config, json_file, indent=4)  # Salvar com formatação (indentação)


    isite, ic = inicializar_rede(Lx, Ly, rc)


    historico_tempo, historico_theta = simular(Lx, Ly, LT, rA, rD, isite, ic, intervalo_snapshot, dirname, modo=modo)




    #subprocess.run(command, text=True)
    print("Simulação completa. Resultados salvos no diretório:", dirname)


if __name__ == "__main__":
    nome_arquivo_config = sys.argv[1] 
    main(nome_arquivo_config)

