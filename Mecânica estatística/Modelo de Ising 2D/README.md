# 🧲 Modelo de Ising 2D com Visualização e Transição de Fase

Simulação do modelo de Ising 2D usando o algoritmo de Metrópolis para diferentes temperaturas, com visualização ao vivo da evolução dos spins e análise da magnetização média por spin em função da temperatura.

## 📁 Estrutura do Projeto

```
.
├── gera_experimentos.py        # Gera 30 diretórios com parâmetros variando a temperatura de 1.0 a 4.0
├── sim.py                      # Código da simulação do modelo de Ising
├── rodatudo.py                 # Executa sim.py dentro de todas as pastas simulacao_*/
├── transicao_fase.py           # Agrega os dados de saída e plota a curva de magnetização média
├── run_all.sh                  # Executa todas as etapas (gera, roda e plota)
├── clean.sh                    # Remove pastas simulacao_* e arquivos gerados
└── README.md                   # Este arquivo
```

---

## 🚀 Como Executar

### 1. Executar todos os passos automaticamente

```bash
chmod +x run_all.sh
./run_all.sh
```

Isso irá:

- Criar 30 pastas com `parametros.json`
- Rodar `sim.py` para cada experimento
- Agregar os resultados e gerar a curva da magnetização média por spin vs. temperatura

### 2. Limpar os arquivos gerados

```bash
chmod +x clean.sh
./clean.sh
```

Remove os diretórios de simulação e arquivos gerados (`*.csv`, `*.png`).

---

## 📈 Saída Esperada

- Arquivo: `grafico_magnetizacao_vs_temperatura.png`
- Arquivo: `curva_magnetizacao_vs_temperatura.csv`
- Vídeo `mp4` com a simulação em cada pasta `simulacao_Tx_xx/`
- CSVs com dados de energia, magnetização e magnetização média

---

## ⚙️ Requisitos

- Python 3
- Bibliotecas:
  - numpy
  - matplotlib
  - json
  - os, re, csv

Instale as dependências com:

```bash
pip install numpy matplotlib
```

---

## 📌 Observação

A linha vertical no gráfico de transição de fase marca a **temperatura crítica** do modelo de Ising 2D:

```
T_c ≈ 2.269
```

---

## 🧪 Licença

Licença MIT


