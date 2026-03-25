# Simulação Numérica da Propagação de Pressão em Meio Unidimensional

## Descrição

Este repositório contém implementações numéricas para a simulação da propagação de pressão em uma mangueira unidimensional utilizando duas abordagens computacionais:

* Método de Diferenças Finitas Explícitas para a equação da onda;
* Método Lattice Boltzmann D1Q3 para transporte unidimensional.

O objetivo é comparar duas metodologias numéricas aplicadas ao mesmo problema físico, analisando estabilidade, comportamento temporal e perfil espacial da solução.

---

## Problema físico modelado

A pressão ( p(x,t) ) é modelada pela equação da onda unidimensional:

\frac{\partial^2 p}{\partial t^2}=\alpha^2\frac{\partial^2 p}{\partial x^2}

onde:

* ( x \in [0,1] )
* ( t \in [0,1] )
* ( \alpha = \sqrt{2} )

---

## Condições do problema

### Condições de contorno

p(0,t)=0.5 \quad ; \quad p(1,t)=1.8

### Condição inicial

p(x,0)=\cos(x)\sin(2\pi x)

### Velocidade inicial

\frac{\partial p}{\partial t}(x,0)=\cos(2\pi x)

---


# 1. Método de Diferenças Finitas

## Formulação numérica

O esquema explícito utilizado é:

w_i^{j+1}=2(1-\lambda^2)w_i^j+\lambda^2(w_{i+1}^j+w_{i-1}^j)-w_i^{j-1}

com:

\lambda=\frac{\alpha k}{h}

onde:

* ( h ): passo espacial
* ( k ): passo temporal

## Critério de estabilidade

O método exige:

\lambda \leq 1

No código:

```python
if lam > 1:
    raise ValueError("Método instável: lambda > 1")
```

## Saída produzida

* valor numérico de ( p(0.5,1) )
* superfície 3D da propagação temporal

---

# 2. Método Lattice Boltzmann D1Q3

## Modelo D1Q3

Velocidades discretas:

c_i = {0,1,-1}

Pesos:

w_i = \left{\frac{2}{3},\frac{1}{6},\frac{1}{6}\right}

## Equação de colisão

f_i(x,t+1)=f_i(x,t)-\frac{1}{\tau}(f_i-f_i^{eq})

## Equilíbrio

f_i^{eq}=w_i p

## Streaming

Distribuição propagada para vizinhos:

* direção positiva → direita
* direção negativa → esquerda

## Saída produzida

* valor numérico de ( p(0.5,1) )
* superfície 3D temporal do campo de pressão


---

# Comparação entre os métodos

| Aspecto                         | Diferenças Finitas | LBM D1Q3                |
| ------------------------------- | ------------------ | ----------------------- |
| Formulação                      | PDE clássica       | Método mesoscópico      |
| Estabilidade                    | CFL                | Relaxação τ             |
| Interpretação física            | Direta             | Distribuições discretas |
| Extensão para dimensões maiores | Moderada           | Natural                 |

---

