import os
import subprocess

diretorio_base = "."

for nome_pasta in os.listdir(diretorio_base):
    caminho_pasta = os.path.join(diretorio_base, nome_pasta)
    if os.path.isdir(caminho_pasta) and os.path.isfile(os.path.join(caminho_pasta, "parametros.json")):
        print(f"Executando simulação em: {caminho_pasta}")
        subprocess.run(["python", "sim.py", caminho_pasta])

