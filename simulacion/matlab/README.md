# Flujo MATLAB

MATLAB ensambla y resuelve el modelo de conducción estacionaria de una placa idealizada. `../00_Proyecto_Principal.mlx` es el único punto de entrada didáctico: sigue el problema desde la ecuación hasta los archivos que consume Python.

## Ejecución

Desde la raíz del repositorio:

```matlab
addpath(fullfile(pwd, 'simulacion', 'matlab'))
checks = validar_simulacion();
assert(checks.passed)
experimentos();
report = generar_reporte_validacion();
```

Para abrir el Live Script académico en MATLAB, abre `../00_Proyecto_Principal.mlx` y selecciona **Run**. También puede ejecutarse en modo batch:

```powershell
matlab -batch "run(fullfile(pwd,'simulacion','00_Proyecto_Principal.m'));"
```

## Responsabilidades

- `../configuracion/modelo_base.json` mantiene los parámetros compartidos.
- `planteamiento_modelo.m` y `configuracion_modelo.m` leen y validan los parámetros.
- `discretizacion_matriz.m` construye la malla, fuentes, matriz dispersa \(A\) y vector \(b\), incorporando Dirichlet.
- `metodo_gauss_seidel.m` implementa Gauss--Seidel lexicográfico con la convención del curso \(A=D-L-U\) y una resolución triangular dispersa.
- `metodo_gradiente_conjugado.m` implementa las recurrencias explícitas de Gradiente Conjugado; no usa `pcg`.
- `comparacion_metodos.m` presenta iteraciones, residuo y temperatura máxima frente a la referencia directa.
- `experimentos.m` ejecuta el caso base, convergencia, refinamiento, sensibilidad, estrés y exportaciones.
- `validar_simulacion.m` comprueba potencias, frontera, simetría, definitud positiva, convergencia, linealidad y referencia directa.
- `generar_reporte_validacion.m` reúne los resultados en un reporte JSON.

El residuo relativo es el criterio común de parada. GS también registra el error relativo entre iteraciones como indicador complementario; este no altera la tolerancia ni las iteraciones reportadas.

## Exportaciones

Los datos se escriben en `../resultados/` como CSV, JSON y MAT. El reporte de validación se guarda en `../pruebas/reporte_validacion.json`. El módulo Python solo lee estas salidas para crear figuras. Los tiempos se reportan por repetición y dependen del equipo y de la carga de la sesión.
