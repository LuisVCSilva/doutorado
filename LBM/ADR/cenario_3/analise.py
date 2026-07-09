import numpy as np
import pandas as pd
from pathlib import Path

print("="*100)
print("ANÁLISE COMPARATIVA COMPLETA - CENÁRIO 2")
print("="*100)

results = []
base = Path("dados")

# Mapeamento completo
solvers = {
    "D1Q2":   base / "d1q2_ruim" / "solution_1d.csv",
    "D1Q3":   base / "d1q3_bom" / "solution_1d.csv",
    "D2Q5":   base / "d2q5_ruim" / "solution.csv",
    "D2Q9":   base / "d2q9_bom" / "solution.csv",
    "FDM_1D": base / "fdm1d_ruim" / "solution.csv",
    "FDM_2D": base / "fdm_adr_2d_results" / "solution.csv",
}

for name, path in solvers.items():
    if path.exists():
        try:
            df = pd.read_csv(path)
            
            # Detectar coluna da solução numérica
            if 'phi_LBM' in df.columns:
                phi = df['phi_LBM'].values
            elif 'phi_fdm' in df.columns:
                phi = df['phi_fdm'].values
            elif 'phi' in df.columns:
                phi = df['phi'].values
            else:
                continue
                
            phi_exact = df['phi_exact'].values
            
            metrics = {
                'Método': name,
                'L2': np.sqrt(np.sum((phi - phi_exact)**2) / np.sum(phi_exact**2)),
                'L∞': np.max(np.abs(phi - phi_exact)),
                'RMSE': np.sqrt(np.mean((phi - phi_exact)**2)),
                'MAE': np.mean(np.abs(phi - phi_exact)),
                'Pontos': len(phi)
            }
            results.append(metrics)
            
            print(f"{name:8} | L2 = {metrics['L2']:.2e} | L∞ = {metrics['L∞']:.2e} | RMSE = {metrics['RMSE']:.2e}")
            
        except Exception as e:
            print(f"Erro ao ler {name}: {e}")
    else:
        print(f"Não encontrado: {name}")

if results:
    df = pd.DataFrame(results)
    df = df.sort_values('L2')
    
    print("\n" + "="*100)
    print("RANKING FINAL (menor L2 = melhor)")
    print("="*100)
    print(df.round(6).to_string(index=False))
    
    df.to_csv("relatorio_completo_cenario2.csv", index=False)
    print(f"\nRelatório salvo em: relatorio_completo_cenario2.csv")
