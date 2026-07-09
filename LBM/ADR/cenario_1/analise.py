import numpy as np
import pandas as pd
from pathlib import Path

print("="*90)
print("ANÁLISE COMPARATIVA - CENÁRIO 1")
print("="*90)

results = []
base = Path("dados")

# Busca automática por arquivos solution
for csv_file in base.rglob("solution*.csv"):
    method = csv_file.parent.name.split('_')[0].upper()
    if "d1q2" in str(csv_file).lower():
        method = "D1Q2"
    elif "d1q3" in str(csv_file).lower():
        method = "D1Q3"
    elif "d2q9" in str(csv_file).lower():
        method = "D2Q9"
    elif "d2q5" in str(csv_file).lower():
        method = "D2Q5"
    elif "fdm" in str(csv_file).lower():
        method = "FDM"

    try:
        df = pd.read_csv(csv_file)
        
        if 'phi_LBM' in df.columns:
            phi = df['phi_LBM'].values
        elif 'phi_fdm' in df.columns:
            phi = df['phi_fdm'].values
        elif 'phi' in df.columns:
            phi = df['phi'].values
        else:
            continue
            
        if 'phi_exact' not in df.columns:
            continue
            
        phi_exact = df['phi_exact'].values
        
        metrics = {
            'method': method,
            'L2': np.sqrt(np.sum((phi - phi_exact)**2) / np.sum(phi_exact**2)),
            'Linf': np.max(np.abs(phi - phi_exact)),
            'RMSE': np.sqrt(np.mean((phi - phi_exact)**2)),
            'MAE': np.mean(np.abs(phi - phi_exact)),
            'points': len(phi),
            'file': str(csv_file)
        }
        results.append(metrics)
        print(f"{method:8} | L2 = {metrics['L2']:.2e} | L∞ = {metrics['Linf']:.2e}")
        
    except:
        pass

if results:
    df_final = pd.DataFrame(results)
    df_final = df_final.sort_values('L2')
    print("\n" + "="*90)
    print("RANKING FINAL")
    print("="*90)
    print(df_final[['method', 'L2', 'Linf', 'RMSE', 'points']].round(6).to_string(index=False))
    
    df_final.to_csv("relatorio_comparativo_cenario1.csv", index=False)
else:
    print("Nenhum arquivo solution*.csv encontrado.")
    print("Execute os solvers primeiro!")
