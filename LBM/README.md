# Lattice Boltzmann Examples

Este repositório contém três implementações simples do método de **Lattice Boltzmann Method** em Python, ilustrando diferentes níveis de complexidade:

* D1Q3 unidimensional
* D2Q9 bidimensional
* D2Q9 multifásico com **Shan–Chen model**

Todos os códigos utilizam:

* colisão BGK
* streaming explícito
* cálculo de grandezas macroscópicas
* exportação de frames PNG

---

## 1. D1Q3 — Difusão unidimensional

Modelo mínimo com três velocidades discretas:

$$
c = {0, +1, -1}
$$

Simula a evolução de uma perturbação de densidade em 1D.

### Física observada

* propagação de ondas de densidade
* relaxação ao equilíbrio
* difusão discreta

### Saída

Frames PNG em:

```bash
lbm_frames/
```

---

## 2. D2Q9 — Difusão bidimensional

Modelo bidimensional com nove velocidades discretas (**D2Q9 lattice**).

A condição inicial consiste em uma perturbação central de densidade.

### Física observada

* expansão isotrópica da perturbação
* relaxação viscosa
* aproximação do regime hidrodinâmico

### Saída

Frames PNG em:

```bash
lbm_d2q9_frames/
```

---

## 3. D2Q9 multifásico — Shan-Chen

Extensão multifásica usando força pseudopotencial:

$$
F(x)=-G,\psi(x)\sum_i w_i \psi(x+c_i)c_i
$$

com:

$$
\psi(\rho)=1-e^{-\rho}
$$

### Física observada

* separação espontânea de fases
* formação de interfaces difusas
* decomposição spinodal

### Saída

Frames PNG em:

```bash
lbm_multiphase/
```

---

## Dependências

```bash
pip install numpy matplotlib
```

---

## Gerar vídeo MP4

Usar **FFmpeg**:

```bash
ffmpeg -framerate 20 -i frame_%04d.png -c:v libx264 -pix_fmt yuv420p output.mp4
```

---

## Observação

Os códigos têm objetivo didático e minimalista. Não incluem:

* condições de contorno complexas
* estabilidade avançada
* conservação rigorosa para casos extremos


