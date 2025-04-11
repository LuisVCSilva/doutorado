import os
import re
import csv
import matplotlib.pyplot as plt

resultados = []


for pasta in os.listdir():
    if os.path.isdir(pasta) and pasta.startswith("simulacao_"):
        arquivos = os.listdir(pasta)
        for arquivo in arquivos:
            if arquivo.endswith("magnetizacao_media.csv"):

                match = re.search(r"T(\d+)_(\d+)", arquivo)
                if match:
                    temp_str = f"{match.group(1)}.{match.group(2)}"
                    temperatura = float(temp_str)

                    caminho_csv = os.path.join(pasta, arquivo)
                    with open(caminho_csv, "r") as f:
                        leitor = csv.reader(f)
                        next(leitor) 
                        linha = next(leitor)
                        magnetizacao_media = float(linha[0])
                        resultados.append((temperatura, magnetizacao_media))

resultados.sort()


with open("curva_magnetizacao_vs_temperatura.csv", "w", newline="") as f:
    escritor = csv.writer(f)
    escritor.writerow(["Temperatura", "MagnetizacaoMedia"])
    escritor.writerows(resultados)


temperaturas, magnetizacoes = zip(*resultados)

plt.figure(figsize=(8, 5))
plt.plot(temperaturas, magnetizacoes, marker='o', linestyle='-', color='blue')
plt.axvline(x=2.269, color='red', linestyle='--', label='Temperatura Crítica $T_c \\approx 2.269$')
plt.title("Magnetização Média por Spin vs Temperatura")
plt.xlabel("Temperatura")
plt.ylabel("Magnetização Média por Spin")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("grafico_magnetizacao_vs_temperatura.png")
plt.show()

