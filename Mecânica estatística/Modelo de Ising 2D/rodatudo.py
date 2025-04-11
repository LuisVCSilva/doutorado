import os
import subprocess

pastas = [p for p in os.listdir() if os.path.isdir(p) and p.startswith("simulacao_")]

for pasta in sorted(pastas):
    caminho_script = os.path.join(pasta, "sim.py")
    if os.path.isfile(caminho_script):
        print(f"Executando: {caminho_script}")
        try:
            subprocess.run(["python", "sim.py"], cwd=pasta, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Erro ao executar {caminho_script}: {e}")
    else:
        print(f"Script não encontrado em {pasta}")
