# THERMOPLATE

THERMOPLATE resuelve la conducción de calor estacionaria en una placa rectangular idealizada con dos fuentes internas de generación, asociadas al modelo con regiones CPU y GPU. El repositorio comparte únicamente la implementación computacional, sus pruebas y las instrucciones para reproducirla. Los datos de simulación y las figuras se generan localmente y no forman parte del repositorio.

## Modelo

Se resuelve la ecuación de Poisson

\[
\nabla^2 T(x,y) + \frac{q'''(x,y)}{k} = 0 \quad \text{en } \Omega,
\qquad T = T_b \quad \text{en } \partial\Omega.
\]

Un stencil de cinco puntos de segundo orden transforma el problema en un sistema disperso \(A\mathbf T=\mathbf b\). Los nodos del contorno de Dirichlet se fijan y sus contribuciones se incorporan en el vector derecho. Las fuentes rectangulares se integran con pesos trapezoidales para conservar la potencia discretizada.

El escenario configurado en `simulacion/configuracion/modelo_base.json` usa una placa de 0.24 × 0.16 m, espesor de 0.002 m, conductividad de 50 W/(m·K) y contorno a 25 °C. Las huellas CPU/GPU miden 0.04 × 0.04 m y reciben 12 W y 18 W. Son supuestos ilustrativos: el modelo no representa la temperatura real de una computadora y no incluye convección, radiación, ventilación, contacto térmico, materiales multicapa ni transitorios.

## Solucionadores

- **Gauss–Seidel:** actualización lexicográfica propia sobre el sistema disperso.
- **Gradiente Conjugado (CG):** recurrencias propias de residuo y dirección conjugada; el código no delega a `pcg`.
- **Referencia directa:** división dispersa de MATLAB, utilizada para medir el error algebraico.

La tolerancia relativa del caso base es \(10^{-8}\). Los experimentos incluyen escenario base, convergencia, cuatro mallas, sensibilidad de cada fuente y pruebas de estrés. El residuo algebraico, la diferencia entre mallas y las simplificaciones físicas se interpretan por separado.

## Requisitos

- MATLAB R2025a para resolver y validar las simulaciones.
- Python 3.10 o posterior, NumPy, Pandas y Matplotlib para generar figuras desde las exportaciones MATLAB.
- Git para clonar y versionar el código.

Instala las dependencias Python desde la raíz del proyecto:

```powershell
python -m pip install -r simulacion/python/requirements.txt
```

No se requiere SciPy para el flujo recomendado: MATLAB exporta CSV y JSON, que son los formatos leídos por las figuras.

## Reproducción en MATLAB R2025a

Abre la carpeta del repositorio en MATLAB o ejecútalo desde una terminal situada en la raíz. Este comando valida el modelo y ejecuta todos los experimentos:

```powershell
matlab -batch "addpath(fullfile(pwd,'simulacion','matlab')); checks=verify_thermoplate_core(); assert(checks.passed); run_thermoplate_experiments(); report=write_thermoplate_validation_report();"
```

En MATLAB, también puedes ejecutar las funciones por separado:

```matlab
addpath(fullfile(pwd, 'simulacion', 'matlab'))
checks = verify_thermoplate_core();
assert(checks.passed)
results = run_thermoplate_experiments();
report = write_thermoplate_validation_report();
```

Los Live Scripts `simulacion/matlab/01_escenario_base.mlx` y `02_experimentos_completos.mlx` documentan la ejecución desde MATLAB. El experimento completo puede tardar varios minutos; en la ejecución de referencia tomó cerca de tres minutos y depende del equipo.

## Reproducción en Python desde VS Code o PyCharm

Abre la raíz del repositorio como carpeta o proyecto, selecciona un intérprete Python 3.10 o posterior y abre la terminal integrada. Instala las dependencias y ejecuta primero MATLAB para crear las exportaciones. Después genera las figuras:

```powershell
python -m pip install -r simulacion/python/requirements.txt
python simulacion/python/visualizar.py --strict --results simulacion/resultados --output simulacion/figuras
```

Python solo lee los resultados exportados; no crea temperaturas ficticias ni vuelve a resolver el modelo térmico. `--strict` devuelve error si faltan datos o metadatos requeridos. Para comprobar el lector y la validación de entradas:

```powershell
python -B -m unittest discover -s simulacion/python/tests -v
```

## Archivos generados

MATLAB crea bajo `simulacion/resultados/` tablas CSV, manifiestos y metadatos JSON, tiempos por repetición y el archivo `baseline.mat` para uso posterior en MATLAB. Los escenarios son `escenario_base`, `convergencia`, `refinamiento`, `sensibilidad` y `estres_numerico`. El reporte de validación se escribe en `simulacion/pruebas/reporte_validacion.json`.

Python lee esas exportaciones y escribe nueve pares de figuras PNG (300 dpi) y PDF vectorial en `simulacion/figuras/`. Tanto esos resultados como las figuras están excluidos de Git; se vuelven a crear al ejecutar los comandos anteriores.

## Resultados de referencia

La ejecución MATLAB R2025a registrada obtuvo para la malla base de 49 × 33 nodos (1,457 incógnitas) un máximo de 77.7592 °C. Gauss–Seidel convergió en 2,145 iteraciones y CG en 127, ambos con residuo relativo inferior a \(10^{-8}\). La matriz fue simétrica, el error de frontera fue cero y las potencias discretizadas fueron 12 W y 18 W. Estos valores son comprobaciones de referencia; al repetir el flujo, consulta los archivos generados para la evidencia de esa ejecución.

## Estructura del repositorio

```text
README.md
simulacion/
  configuracion/   Parámetros de simulación compartidos
  matlab/          Ensamblaje, solucionadores, experimentos y validación
  python/          Visualización, dependencias y pruebas automatizadas
```

Los directorios `simulacion/resultados/`, `simulacion/figuras/` y `simulacion/pruebas/` aparecen después de ejecutar MATLAB y Python; su contenido generado no se versiona. La carpeta académica `Latex/` se conserva localmente y está excluida del repositorio.
