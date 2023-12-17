from tqdm import tqdm 
import numpy as np
import matplotlib.pyplot as plt

const_boltzmann = 1.3806488e-23 
hbar = 6.626070040e-34  
temp_rot = 15.2 #Para HCl = 15.2 Outros valores em P. Atkins and J. de Paula "Physical Chemistry", 9th edition (W.H. Freeman 2010), Table 13.2, Data section in appendix

def energia(l):
    return const_boltzmann * l * (l + 1) * temp_rot

def integracao_monteCarlo(temperatura, num_steps=5000, num_trials=10000):
    energias = np.zeros(num_trials)
    quantum_l = np.zeros(num_trials)
    
    for c2 in tqdm(range(num_trials), desc=f'temperatura: {temperatura} K', unit='trial'):
        e_sum = 0
        e_num = 0
        l = 10 

        for c1 in range(num_steps):
            lz = l + 1 if np.random.rand() < 0.5 else l - 1 

            if lz < 0:
                p = 0
            else:
                p = (2.0 * lz + 1.0) / (2.0 * l + 1.0) * np.exp(-(energia(lz) - energia(l)) / (const_boltzmann * temperatura))

            if np.random.rand() < p:
                l = lz

            e_sum += energia(l)
            e_num += 1

        energias[c2] = e_sum / e_num
        quantum_l[c2] = l

    return energias, quantum_l


def main():
    temperaturas = np.arange(1, 61)
    energias_medias = []

    for temperatura in temperaturas:
        energias, _ = integracao_monteCarlo(temperatura)
        energias_medias.append(np.mean(energias))

    plt.plot(temperaturas, energias_medias, marker='o')
    plt.xlabel('Temperatura')
    plt.ylabel('Energia média')
    plt.show()

if __name__ == "__main__":
    main()
