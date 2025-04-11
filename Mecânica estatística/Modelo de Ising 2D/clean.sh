#!/bin/bash

echo "Limpando diretórios de simulação..."
rm -rf simulacao_*/

echo "Removendo arquivos de saída gerais..."
rm -f curva_magnetizacao_vs_temperatura.csv
rm -f grafico_magnetizacao_vs_temperatura.png

echo "Limpeza concluída!"
