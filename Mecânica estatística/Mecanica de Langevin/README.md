# Equação de Langevin e Simulação com Leapfrog

Este repositório contém códigos Python para simular o movimento de uma partícula sob a ação de forças determinísticas e estocásticas, descritas pela **equação de Langevin**. A simulação é feita utilizando o método **Leapfrog com meio passo**, que possui vantagens numéricas para sistemas com forças conservativas e dissipativas.

---

## 📘 Equação de Langevin

A **equação de Langevin** descreve a evolução do momento linear de uma partícula sujeita a três tipos de forças:

$$
\frac{d\vec{p}}{dt} = -\beta \vec{p} + \vec{F}_{\text{ex}}(\vec{r}) + \vec{f}(t)
$$

Onde:
- \( \vec{p}(t) \): momento linear da partícula.
- \( \beta \): coeficiente de dissipação (atrito viscoso).
- \( \vec{F}_{\text{ex}}(\vec{r}) \): força externa determinística (ex: força de mola).
- \( \vec{f}(t) \): força aleatória (ruído térmico).

---

## 🎲 Propriedades Estatísticas do Ruído Aleatório \( \vec{f}(t) \)

O ruído é modelado como **ruído branco gaussiano**, com as seguintes propriedades:

- **Média nula**:
  $$
  \langle f_x(t) \rangle = \langle f_y(t) \rangle = \langle f_z(t) \rangle = 0
  $$

- **Sem correlação temporal**:
  $$
  \langle f_i(t) f_j(t') \rangle = F_0^2 \delta_{ij} \delta(t - t')
  $$

- **Variância finita** em cada direção:
  $$
  \langle f_i^2(t) \rangle = F_0^2
  $$

- **Amplitude total das flutuações**:
  $$
  \langle |\vec{f}(t)|^2 \rangle = 3 F_0^2
  $$

---

## 💻 Código: Leapfrog com Meio Passo

O código implementa o método Leapfrog, onde as **velocidades são calculadas em meio passo**:

1. **Inicialização** com meio passo:
   $$
   v_{i+\frac{1}{2}} = v_i + \frac{dt}{2} \cdot a_i
   $$

2. **Evolução no tempo**:
   - Atualiza posição: \( x_{i+1} = x_i + dt \cdot v_{i+\frac{1}{2}} \)
   - Calcula nova aceleração: inclui dissipação, força restauradora e ruído
   - Atualiza velocidade de meio passo: \( v_{i+\frac{3}{2}} = v_{i+\frac{1}{2}} + dt \cdot a_{i+1} \)

---

## 📎 Requisitos

- Python 3
- NumPy
- Matplotlib

Instale com:

```bash
pip install numpy matplotlib
