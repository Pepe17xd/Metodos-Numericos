# Modelación numérica de la distribución térmica en una placa electrónica idealizada mediante Diferencias Finitas, Gauss–Seidel y Gradiente Conjugado

Proyecto del curso CC2104, Métodos Numéricos (UTEC). Se calcula el campo estacionario de temperatura en una placa rectangular idealizada y se comparan Gauss–Seidel y Gradiente Conjugado sobre el sistema obtenido por diferencias finitas.

## Modelo

Se resuelve

\[
\nabla^2 T + \frac{q'''(x,y)}{k}=0\quad\text{en }\Omega,
\qquad T=T_b\quad\text{en }\partial\Omega.
\]

La CPU y la GPU representan dos regiones rectangulares de generación uniforme. Sus potencias de 12 W y 18 W, así como las propiedades de la placa, son supuestos de simulación. No representan mediciones de una computadora. El artículo científico base es Santos y Gomes dos Anjos (2022), *Finite Difference Method Applied in Two-Dimensional Heat Conduction Problem in the Permanent Regime in Rectangular Coordinates*; sus ecuaciones se usan como referencia matemática, mientras que la geometría, los parámetros, los experimentos y la comparación de solucionadores corresponden a esta adaptación.

Las diferencias finitas ensamblan el sistema disperso \(A\mathbf{T}=\mathbf{b}\). Gauss–Seidel realiza barridos secuenciales; Gradiente Conjugado usa direcciones conjugadas para una matriz simétrica definida positiva. Una solución directa dispersa de MATLAB sirve como referencia de verificación. La parada de ambos métodos iterativos usa el mismo residuo relativo \(\|b-Ax\|_2/\|b\|_2\).

## Organización

```text
simulacion/
  00_Proyecto_Principal.mlx  Laboratorio guiado: físico -> ecuación -> A*x=b -> GS -> CG -> resultados
  00_Proyecto_Principal.m    Misma secuencia en formato de script para ejecución por lotes
  matlab/                    Funciones breves de modelo, ensamblaje, métodos, experimentos y validación
  python/graficas.py         Figuras a partir de exportaciones MATLAB; no resuelve el modelo
  resultados/                CSV, JSON y MAT exportados por MATLAB
  figuras/                   Figuras PNG/PDF generadas desde esas exportaciones
  pruebas/                   Reporte JSON de validación de MATLAB
```

Abra primero `simulacion/00_Proyecto_Principal.mlx`. El laboratorio sigue el orden del curso: problema físico, ecuación del calor, diferencias finitas, ejemplo de \(A\mathbf{x}=\mathbf{b}\), Gauss--Seidel, Gradiente Conjugado, comparación y exportación para Python. Las funciones de `matlab/` mantienen cada responsabilidad visible y usan matrices dispersas; no forman inversas ni llaman a `pcg`.

## Requisitos

- MATLAB R2025a Update 1 o posterior para resolver y validar el modelo. La validación incluida se ejecutó con R2025a Update 1.
- Python 3.10 o posterior, NumPy, Pandas y Matplotlib para visualizar las exportaciones.
- Git para descargar el proyecto.

Instala las dependencias de visualización desde la raíz:

```powershell
python -m pip install -r simulacion/python/requirements.txt
```

## Reproducir

Abre el Live Script principal en MATLAB y ejecútalo de principio a fin, o usa el modo batch:

```powershell
matlab -batch "run(fullfile(pwd,'simulacion','00_Proyecto_Principal.m'));"
```

El flujo ejecuta el escenario base, refinamiento de malla, sensibilidad a las potencias, pruebas de tolerancia y exportación. Para correr solo la validación central:

```powershell
matlab -batch "addpath(fullfile(pwd,'simulacion','matlab')); checks=validar_simulacion(); disp(checks); assert(checks.passed);"
```

Para repetir únicamente los experimentos adicionales y actualizar `simulacion/resultados/`:

```powershell
matlab -batch "addpath(fullfile(pwd,'simulacion','matlab')); experimentos();"
```

Después de MATLAB, genera las gráficas a partir de los archivos exportados:

```powershell
python simulacion/python/graficas.py --strict --results simulacion/resultados --output simulacion/figuras
```

Las figuras cubren el mapa térmico 2D, la superficie 3D, estados iterativos, residuos, tiempos por refinamiento, sensibilidad de temperatura y pruebas de estrés. Python no resuelve el sistema ni crea datos térmicos. Las pruebas del lector se ejecutan con:

```powershell
python -B -m unittest discover -s simulacion/python/tests -v
```

## Resultados y alcance

Los CSV, JSON y MAT de `simulacion/resultados/` contienen campos de temperatura, historiales, estados representativos, resúmenes y tiempos por repetición; `simulacion/pruebas/reporte_validacion.json` registra las comprobaciones. Las figuras de `simulacion/figuras/` presentan esos datos exportados. Las salidas versionadas son una instantánea reproducible y se pueden regenerar con MATLAB y Python.

En la ejecución de referencia registrada, la malla de 49 × 33 nodos produjo \(T_{\max}=77.7592\ ^\circ\mathrm{C}\); Gauss–Seidel alcanzó residuo menor que \(10^{-8}\) en 2145 iteraciones y Gradiente Conjugado en 127. Los tiempos cambian con el equipo y la carga del sistema. El modelo es estacionario, bidimensional, homogéneo e isotrópico, con conductividad constante y temperatura prescrita en todo el borde. No incorpora convección, radiación, ventiladores, contacto, transitorios ni capas de una placa real.
