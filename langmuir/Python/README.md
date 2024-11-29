python3 -m venv isoterma # cria ambiente virtual
source isoterma/bin/activate # ativa ambiente virtual

pip install -r requirements.txt # instala dependências
python3 isoterma.py cenarios/favoravel.json # roda KMC com isoterma favorável


rm -rf $(ls | grep sim_) #Limpar simulações
