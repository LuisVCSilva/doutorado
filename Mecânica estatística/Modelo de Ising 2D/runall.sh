#!/bin/bash

echo "Rodando gera_experimentos.py..."
python3 gera_experimentos.py

echo "Rodando rodatudo.py..."
python3 rodatudo.py

echo "Rodando transicao_fase.py..."
python3 transicao_fase.py

echo "Processo completo!"
