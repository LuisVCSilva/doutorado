from numpy import *
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import ipywidgets as widgets
from IPython.display import display

# constantes da membrana
C_M = 1.0           # capacitancia da membrana por unidade de area (F/m^2)
V_Na = -115         # potencial nerst do sodio (mV)
V_K = +12           # potencial nerst do potassio (mV)
V_l = -10.613       # potencial de vazamento (mv)
g_Na = 120          # condutancia do sodio (S/m^2)
g_K = 36            # condutancia do potassio (S/m^2)
g_l = 0.3           # condutancia do vazamento (S/m^2)

# constantes de taxa dependentes de voltagem 
def alpha_n(V):
    return 0.01 * (V + 10) / (exp((V + 10) / 10) - 1)

def beta_n(V):
    return 0.125 * exp(V / 80)

def alpha_m(V):
    return 0.1 * (V + 25)  / (exp((V + 25) / 10) - 1)

def beta_m(V):
    return 4 * exp(V / 18)

def alpha_h(V):
    return 0.07 * exp(V / 20)

def beta_h(V):
    return 1 / (exp((V + 30) / 10) + 1)

def I(t):
    return 0

def Modelo_HH(t, y):
    V, n, m, h = y
    saida = [0] * 4
    saida[0] = ( I(t) - g_K * n**4 * (V - V_K) - g_Na * m**3 * h * (V - V_Na) - g_l * (V - V_l) ) / C_M
    saida[1] = alpha_n(V) * (1 - n) - beta_n(V) * n
    saida[2] = alpha_m(V) * (1 - m) - beta_m(V) * m
    saida[3] = alpha_h(V) * (1 - h) - beta_h(V) * h
    return saida

    # Corrente da membrana


n_0 = alpha_n(0) / (alpha_n(0) + beta_n(0))
m_0 = alpha_m(0) / (alpha_m(0) + beta_m(0))
h_0 = alpha_h(0) / (alpha_h(0) + beta_h(0))
V_0 = -90
y_0 = [V_0, n_0, m_0, h_0]
t_span = (0,4)
sol = solve_ivp(Modelo_HH, t_span, y_0, method='RK45', t_eval=linspace(*t_span, 100))

plt.figure()
plt.plot(sol.t, sol.y[0], label='V(t)')
plt.xlabel('Tempo (ms)')
plt.ylabel('Potencial da membrana (mV)')
plt.legend()
plt.show()

plt.figure()
plt.plot(sol.t, sol.y[1], label='n(t)')
plt.plot(sol.t, sol.y[2], label='m(t)')
plt.plot(sol.t, sol.y[3], label='h(t)')
plt.xlabel('Tempo (ms)')
plt.ylabel('Condutância')
plt.legend()
plt.show()
