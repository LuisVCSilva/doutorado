import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime


start_date = datetime(2023, 10, 1)
end_date = datetime(2027, 10, 1)


tasks = [
    ("Ano 1", datetime(2023, 10, 1), datetime(2024, 10, 1)),
    ("Estudo do estado da arte e replicação de resultados", datetime(2023, 10, 1), datetime(2024, 1, 1)),
    ("Escrita do projeto de doutorado", datetime(2024, 1, 1), datetime(2024, 4, 1)),
    ("Investigação de aspectos teóricos do problema-alvo", datetime(2024, 4, 1), datetime(2024, 4, 15)),
    ("Estudo de técnicas teórico-computacionais para resolução de casos clássicos", datetime(2024, 4, 15), datetime(2024, 6, 1)),
    ("Criação de casos para validação", datetime(2024, 6, 1), datetime(2024, 7, 1)),
    ("Investigações experimentais", datetime(2024, 7, 1), datetime(2024, 8, 1)),
    ("Artigo para avaliação de métodos teóricos-computacionais relacionado ao problema-alvo", datetime(2024, 8, 1), datetime(2024, 9, 1)),
    ("Artigo de revisão", datetime(2024, 9, 1), datetime(2024, 10, 1)),
    ("Ano 2", datetime(2024, 10, 1), datetime(2025, 10, 1)),
    ("Organização de dados dos experimentos práticos", datetime(2024, 10, 1), datetime(2025, 1, 1)),
    ("Escrita de artigo experimental", datetime(2025, 1, 1), datetime(2025, 2, 1)),
    ("Organização de resultados e balanço do projeto", datetime(2025, 2, 1), datetime(2025, 3, 1)),
    ("Escrita de artigo principal", datetime(2025, 3, 1), datetime(2025, 4, 1)),
    ("Submissão de resultados principais", datetime(2025, 4, 1), datetime(2025, 6, 1)),
    ("Ano 3", datetime(2025, 10, 1), datetime(2026, 10, 1)),
    ("Escrita da tese (caso não escandinavo)", datetime(2026, 4, 1), datetime(2027, 4, 1)),
    ("Discussões com pares externos", datetime(2027, 4, 1), datetime(2027, 5, 1)),
    ("Revisões internas", datetime(2027, 5, 1), datetime(2027, 6, 1)),
    ("Defesa da tese de doutorado final", datetime(2027, 6, 1), datetime(2027, 6, 15)),
]


fig, ax = plt.subplots()


for i, (task, start, end) in enumerate(tasks):
    ax.barh(y=i, width=end-start, left=start, height=0.5, color='b', alpha=0.8)
    ax.text(start + (end - start) / 2, i, "", ha='center', va='center', color='white')


ax.set_xlim(start_date, end_date)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax.xaxis.set_major_locator(mdates.YearLocator())


ax.set_xlabel('Tempo')
ax.set_yticks(range(len(tasks)))
ax.set_yticklabels([task for task, _, _ in tasks])
ax.set_title('Planejamento de tarefas do doutorado')


plt.xticks(rotation=45)


plt.show()

