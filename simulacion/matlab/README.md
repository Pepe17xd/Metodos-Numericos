# Flujo MATLAB

MATLAB ensambla y resuelve el modelo de conducción estacionaria de una placa idealizada. `../00_Proyecto_Principal.mlx` es el Live Script didáctico principal; sigue el problema desde la ecuación hasta los archivos que consume Python. `02_experimentos_completos.mlx` ejecuta por separado el conjunto de experimentos.

## Ejecución

Desde la raíz del repositorio:

```matlab
addpath(fullfile(pwd, 'simulacion', 'matlab'))
checks = verify_thermoplate_core();
assert(checks.passed)
ejecutar_experimentos();
report = write_thermoplate_validation_report();
```

Para abrir el Live Script académico en MATLAB, abre `../00_Proyecto_Principal.mlx` y selecciona **Run**. También puede ejecutarse en modo batch:

```powershell
matlab -batch "run(fullfile(pwd,'simulacion','00_Proyecto_Principal.m'));"
```

## Responsabilidades

- `../configuracion/modelo_base.json` mantiene los parámetros compartidos.
- `thermoplate_config.m`, `load_thermoplate_config.m` validan la configuración.
- `02_discretizacion_matriz.m` construye la malla, fuentes, matriz dispersa \(A\) y vector \(b\), incorporando Dirichlet.
- `03_metodo_gauss_seidel.m` implementa Gauss--Seidel lexicográfico con la convención del curso \(A=D-L-U\) y una resolución triangular dispersa.
- `04_metodo_gradiente_conjugado.m` implementa las recurrencias explícitas de Gradiente Conjugado; no usa `pcg`.
- `06_experimentos.m` ejecuta el caso base, convergencia, refinamiento, sensibilidad, estrés y exportaciones.
- `verify_thermoplate_core.m` comprueba potencias, frontera, simetría, definitud positiva, convergencia, linealidad y referencia directa.
- `write_thermoplate_validation_report.m` reúne los resultados en un reporte JSON.

El residuo relativo es el criterio común de parada. GS también registra el error relativo entre iteraciones como indicador complementario; este no altera la tolerancia ni las iteraciones reportadas. Los cuatro adaptadores con nombre histórico solo preservan compatibilidad con el Live Script previamente guardado; el flujo documentado utiliza los archivos numerados.

## Exportaciones

Los datos se escriben en `../resultados/` como CSV, JSON y MAT. El reporte de validación se guarda en `../pruebas/reporte_validacion.json`. El módulo Python solo lee estas salidas para crear figuras. Los tiempos se reportan por repetición y dependen del equipo y de la carga de la sesión.
