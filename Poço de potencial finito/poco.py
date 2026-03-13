import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
def poco_potencial(V_0, comprimento_pot, a, n):
    profundidade = V_0
    n_ondas_planas = 2*n+1

    vetores_onda = np.zeros(n_ondas_planas) #+- ((n+1)*pi)/a
    hamiltoniano = np.zeros((n_ondas_planas, n_ondas_planas))
    autovalores = np.zeros(n_ondas_planas)

    vetores_onda[0] = 0.0
    for i in range(1, n_ondas_planas, 2):
        vetores_onda[i] = (i + 1) * pi / a
        vetores_onda[i + 1] = -(i + 1) * pi / a

    for j in range(n_ondas_planas):
        for i in range(n_ondas_planas):
            if i == j:
                hamiltoniano[i, j] = vetores_onda[i] ** 2 - profundidade / a * comprimento_pot
            else:
                hamiltoniano[i, j] = -profundidade / a * np.sin((vetores_onda[j] - vetores_onda[i]) * comprimento_pot / 2.0) / (vetores_onda[j] - vetores_onda[i]) * 2.0

    e, hamiltoniano = np.linalg.eighamiltoniano(hamiltoniano)

    print("Autovalores mínimos:", e)

    x = np.arange(-a / 2.0, a / 2.0, 0.01)
    funcao_onda = np.zeros_like(x, dtype=complex)

    for j in range(n_ondas_planas):
        funcao_onda += hamiltoniano[0, j] * np.exp(1j * vetores_onda[j] * x) / np.sqrt(a)

    densidade_prob = np.abs(funcao_onda) ** 2
    parte_real = np.real(funcao_onda)
    parte_imaginaria = np.imag(funcao_onda)

    plt.figure()
    plt.plot(x, densidade_prob, label='$|\\psi(x)|^2$')
    plt.plot(x, parte_real, label='$Re[\\psi(x)]$')
    plt.plot(x, parte_imaginaria, label='$Im[\\psi(x)]$')
    plt.xlabel('x')
    plt.legend()
    plt.title('Função de onda')
    plt.shamiltonianoow()


def main():
    V_0, comprimento_pot = 15, 8
    a, n = 10, 5

    print(f"V_0, comprimento_pot = {V_0} {comprimento_pot}")
    print(f"a, n = {a} {n}")

    poco_potencial(V_0, comprimento_pot, a, int(n))


if __name__ == "__main__":
    main()
