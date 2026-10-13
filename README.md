
# Modelación numérica de la distribución térmica en una placa electrónica idealizada mediante Diferencias Finitas, Gauss–Seidel y Gradiente Conjugado

**Universidad de Ingeniería y Tecnología (UTEC)**  
**Curso:** CC2104 – Métodos Numéricos  
**Herramientas:** MATLAB y Python

## Descripción del proyecto

Este proyecto estudia la distribución de temperatura en una placa electrónica idealizada, considerando dos fuentes internas de generación de calor que representan una CPU y una GPU.

El objetivo es resolver numéricamente la ecuación de conducción de calor bidimensional en régimen estacionario y comparar el desempeño de dos métodos iterativos:

- **Gauss–Seidel:** obtiene aproximaciones sucesivas mediante actualizaciones secuenciales.
- **Gradiente Conjugado:** utiliza residuos y direcciones conjugadas para resolver sistemas simétricos definidos positivos.

Para construir el sistema de ecuaciones se utiliza el **Método de Diferencias Finitas (MDF)**.

MATLAB realiza los cálculos numéricos y Python genera las gráficas a partir de los resultados exportados.

> **Nota:** las temperaturas y potencias corresponden a un modelo idealizado. No representan mediciones de una laptop real.

---

## Inicio rápido: descargar y ejecutar

**No necesitas saber Git para utilizar este proyecto.**

### Opción 1. Descargar sin Git

1. Ingresa al [repositorio de GitHub](https://github.com/Pepe17xd/Metodos-Numericos).
2. Presiona el botón verde **Code**.
3. Selecciona **Download ZIP**.
4. Extrae el archivo descargado.
5. Abre MATLAB.
6. Busca y abre el archivo:

   `simulacion/00_Proyecto_Principal.mlx`

7. Presiona **Run** para ejecutar el laboratorio.

El proyecto utiliza rutas relativas para localizar sus archivos, por lo que no requiere configurar manualmente rutas absolutas.

**Importante:** la ejecución completa incluye pruebas con mallas grandes y puede tardar varios minutos.

### Opción 2. Descargar usando Git

Abre una terminal y ejecuta:

```bash
git clone https://github.com/Pepe17xd/Metodos-Numericos.git
cd Metodos-Numericos
```

Después abre `simulacion/00_Proyecto_Principal.mlx` en MATLAB.

---

## Requisitos

### Para ejecutar la simulación

- MATLAB (validado con **R2025a Update 1**).
- Repositorio descargado y descomprimido.

### Para generar las gráficas

- Python 3.10 o posterior.
- NumPy.
- Pandas.
- Matplotlib.

Instala las dependencias desde la carpeta principal:

```bash
python -m pip install -r simulacion/python/requirements.txt
```

**Git es opcional:** también puedes descargar el proyecto mediante ZIP.

---

## Organización del proyecto

```text
Metodos-Numericos/
│
├── README.md
│
└── simulacion/
    │
    ├── 00_Proyecto_Principal.mlx
    ├── 00_Proyecto_Principal.m
    │
    ├── matlab/
    │   ├── planteamiento_modelo.m
    │   ├── discretizacion_matriz.m
    │   ├── metodo_gauss_seidel.m
    │   ├── metodo_gradiente_conjugado.m
    │   ├── comparacion_metodos.m
    │   └── experimentos.m
    │
    ├── python/
    │   ├── graficas.py
    │   ├── requirements.txt
    │   └── tests/
    │
    ├── resultados/
    ├── figuras/
    └── pruebas/
```

### ¿Qué hace cada archivo principal?

| Archivo | Función |
|---|---|
| `00_Proyecto_Principal.mlx` | Laboratorio guiado con explicaciones y ejecución |
| `00_Proyecto_Principal.m` | Versión ejecutable mediante comandos |
| `planteamiento_modelo.m` | Define las características físicas de la placa |
| `discretizacion_matriz.m` | Construye el sistema de ecuaciones mediante MDF |
| `metodo_gauss_seidel.m` | Implementa Gauss–Seidel |
| `metodo_gradiente_conjugado.m` | Implementa Gradiente Conjugado |
| `comparacion_metodos.m` | Compara los resultados de ambos solucionadores |
| `experimentos.m` | Ejecuta refinamiento, sensibilidad y pruebas adicionales |
| `graficas.py` | Genera figuras con los resultados de MATLAB |

---

## Metodología numérica

El laboratorio sigue este procedimiento:

**Problema físico → Ecuación del calor → Diferencias Finitas → Sistema lineal → Gauss–Seidel y Gradiente Conjugado → Comparación → Gráficas.**

### 1. Modelo de conducción térmica

La placa se representa mediante un dominio rectangular bidimensional, homogéneo e isotrópico, con conductividad térmica constante.

La ecuación gobernante es:

$$
\frac{\partial^2 T}{\partial x^2}
+
\frac{\partial^2 T}{\partial y^2}
+
\frac{q'''(x,y)}{k}
=0
$$

Donde:

- $T$: temperatura de la placa.
- $q'''(x,y)$: generación volumétrica de calor.
- $k$: conductividad térmica del material.

Se utiliza la condición de frontera:

$$
T=T_b \qquad \text{en } \partial\Omega
$$

Las fuentes internas representan:

- **CPU:** 12 W.
- **GPU:** 18 W.

Estos valores son parámetros asumidos para la simulación.

### 2. Método de Diferencias Finitas

El dominio se divide en una malla de nodos.

Para una malla uniforme con $\Delta x=\Delta y=h$, la temperatura de cada nodo interior se aproxima mediante:

$$
T_{i,j}=
\frac{
T_{i+1,j}
+
T_{i-1,j}
+
T_{i,j+1}
+
T_{i,j-1}
+
\frac{q'''_{i,j}h^2}{k}
}{4}
$$

Al aplicar esta expresión a los nodos interiores se obtiene el sistema:

$$
A\mathbf{x}=\mathbf{b}
$$

La matriz $A$ se construye de forma dispersa para evitar almacenar elementos nulos innecesarios.

### 3. Método de Gauss–Seidel

Se utiliza la descomposición enseñada en el laboratorio del curso:

$$
A=D-L-U
$$

Donde $D$ es la diagonal, mientras que $L$ y $U$ representan los negativos de las partes triangulares inferior y superior estrictas.

La iteración es:

$$
(D-L)x^{(k+1)}=Ux^{(k)}+b
$$

En MATLAB se emplea la resolución triangular, evitando calcular inversas explícitas.

### 4. Método de Gradiente Conjugado

Para la matriz simétrica definida positiva obtenida mediante la discretización, el sistema lineal puede relacionarse con la minimización de:

$$
f(x)=\frac{1}{2}x^TAx-b^Tx
$$

Su gradiente es:

$$
\nabla f(x)=Ax-b
$$

Gradiente Conjugado construye direcciones de búsqueda conjugadas utilizando el residuo y las aproximaciones anteriores.

Se implementa explícitamente sin utilizar la función predefinida `pcg`.

### 5. Criterio de convergencia

Ambos métodos se comparan mediante el residuo relativo:

$$
r_k=
\frac{\|b-Ax^{(k)}\|_2}{\|b\|_2}
$$

La condición de parada es:

$$
r_k \leq 10^{-8}
$$

El método de Gauss–Seidel también registra el error relativo entre aproximaciones sucesivas como indicador complementario.

---

## Ejecución desde la terminal

Los siguientes comandos deben ejecutarse desde la raíz del repositorio.

### Ejecutar el proyecto completo

```powershell
matlab -batch "run(fullfile(pwd,'simulacion','00_Proyecto_Principal.m'));"
```

### Ejecutar la validación central

```powershell
matlab -batch "addpath(fullfile(pwd,'simulacion','matlab')); checks=validar_simulacion(); disp(checks); assert(checks.passed);"
```

### Ejecutar experimentos adicionales

```powershell
matlab -batch "addpath(fullfile(pwd,'simulacion','matlab')); experimentos();"
```

Los experimentos incluyen:

- Refinamiento espacial de la malla.
- Sensibilidad a las potencias de CPU y GPU.
- Comparaciones de convergencia.
- Medición del costo computacional.
- Pruebas de tolerancia y robustez numérica.

---

## Visualizaciones en Python

Una vez ejecutada la simulación en MATLAB, se pueden generar las gráficas con:

```powershell
python simulacion/python/graficas.py --strict --results simulacion/resultados --output simulacion/figuras
```

El programa genera visualizaciones de:

1. Geometría de la placa y fuentes térmicas.
2. Campo de temperatura bidimensional.
3. Superficie de temperatura tridimensional.
4. Estados de las aproximaciones iterativas.
5. Convergencia de Gauss–Seidel y Gradiente Conjugado.
6. Tiempo computacional según el tamaño de malla.
7. Sensibilidad de la temperatura máxima.
8. Pruebas de estrés y tolerancia numérica.
9. Visualizaciones adicionales de comparación.

Las figuras se almacenan en `simulacion/figuras/`.

**Python no resuelve nuevamente el modelo matemático.** Únicamente lee los resultados calculados y exportados por MATLAB.

### Ejecutar pruebas de Python

```powershell
python -B -m unittest discover -s simulacion/python/tests -v
```

---

## Resultados del escenario base

Se utiliza una placa de $0.24 \times 0.16$ m y una malla de **49 × 33 nodos**, equivalente a **1457 incógnitas interiores**.

| Indicador | Resultado |
|---|---:|
| Temperatura máxima | 77.759239 °C |
| Potencia CPU | 12 W |
| Potencia GPU | 18 W |
| Iteraciones Gauss–Seidel | 2145 |
| Iteraciones Gradiente Conjugado | 127 |
| Residuo final Gauss–Seidel | $9.9420\times10^{-9}$ |
| Residuo final Gradiente Conjugado | $9.2837\times10^{-9}$ |
| Tolerancia | $10^{-8}$ |

Los dos métodos iterativos alcanzaron la tolerancia establecida.

En el escenario base, Gradiente Conjugado necesitó menos iteraciones que Gauss–Seidel.

Los tiempos computacionales dependen del equipo y las condiciones de ejecución.

---

## Alcances y limitaciones

El modelo considera:

- Conducción bidimensional estacionaria.
- Material homogéneo e isotrópico.
- Conductividad térmica constante.
- Generación localizada de calor en CPU y GPU.
- Condiciones de temperatura prescrita en el contorno.

No incorpora convección, radiación, ventiladores, disipadores, resistencias de contacto, estructuras multicapa ni comportamiento transitorio.

Por tanto, los resultados describen el comportamiento del modelo idealizado y no predicen directamente las temperaturas de una laptop real.

---

## Artículo científico de referencia

Gonçalves de Brito dos Santos, V., & Gomes dos Anjos, P. T. (2022). Finite difference method applied in two-dimensional heat conduction problem in the permanent regime in rectangular coordinates. *Advances in Pure Mathematics, 12*, 505–518.

**DOI:** https://doi.org/10.4236/apm.2022.129038

El artículo constituye la referencia matemática inicial. La placa electrónica idealizada, las fuentes CPU/GPU y la comparación de solucionadores corresponden a la adaptación desarrollada en este proyecto.
