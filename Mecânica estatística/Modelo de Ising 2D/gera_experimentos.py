import os
import json
import shutil
import numpy as np

parametros_base = {
    "interacao_J": 1.0,
    "campo_H": 0.0,
    "tam_x": 50,
    "tam_y": 50,
    "passos_mc": 500,
    "intervalo_amostra": 1
}

temperaturas = np.linspace(1.0, 4.0, 30)

nome_script = "sim.py"

if not os.path.isfile(nome_script):
    raise FileNotFoundError(f"O arquivo '{nome_script}' não foi encontrado no diretório atual.")

for i, T in enumerate(temperaturas, start=1):
    nome_pasta = f"simulacao_{i:03d}"
    os.makedirs(nome_pasta, exist_ok=True)
    
    parametros = parametros_base.copy()
    parametros["temperatura"] = round(float(T), 5)
    
    caminho_json = os.path.join(nome_pasta, "parametros.json")
    with open(caminho_json, "w") as f:
        json.dump(parametros, f, indent=4)

    shutil.copy(nome_script, nome_pasta)

print("Experimentos gerados com sucesso!")

