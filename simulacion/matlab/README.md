# Simulación MATLAB

Este directorio implementa el modelo estacionario de conducción, los solucionadores explícitos, los experimentos y la exportación de resultados. Requiere MATLAB R2025a; no requiere Optimization Toolbox para los dos métodos iterativos propios.

## Ejecución

Desde la raíz del proyecto:

```matlab
addpath(fullfile(pwd, 'simulacion', 'matlab'))
checks = verify_thermoplate_core();
disp(checks)
run_thermoplate_experiments();
report = write_thermoplate_validation_report();
```

También puede ejecutarse en modo batch desde PowerShell:

```powershell
matlab -batch "addpath(fullfile(pwd,'simulacion','matlab')); checks=verify_thermoplate_core(); disp(checks); run_thermoplate_experiments(); report=write_thermoplate_validation_report();"
```

Los Live Scripts `01_escenario_base.mlx` y `02_experimentos_completos.mlx` presentan la ejecución documentada paso a paso. Las funciones `.m` son la implementación compartida por los scripts y las pruebas.

## Organización

- `thermoplate_config.m`, `load_thermoplate_config.m`: parámetros y validación de configuración.
- `build_thermoplate_model.m`: malla, cuadratura de fuentes, ensamble disperso de \(A\) y \(b\), y grados de libertad interiores.
- `solve_thermoplate_gs.m`, `solve_thermoplate_cg.m`: Gauss–Seidel lexicográfico y Gradiente Conjugado explícitos.
- `thermoplate_field.m`, `thermoplate_metrics.m`: reconstrucción del campo completo y métricas.
- `run_thermoplate_experiments.m`: escenario base, refinamiento, sensibilidad, estrés y exportaciones.
- `verify_thermoplate_core.m`: controles de simetría, positividad, frontera, potencia, linealidad, residuos y fallo esperado por máximo de iteraciones.
- `write_thermoplate_validation_report.m`: reporte JSON consolidado.

## Algoritmos y validación

El sistema usa un stencil de cinco puntos para \(-\nabla^2 T=q'''/k\). Las temperaturas de Dirichlet se fijan en el borde y se trasladan al vector derecho. Las fuentes rectangulares usan pesos trapezoidales en los bordes para que la integral discreta corresponda a la potencia configurada.

Gauss–Seidel actualiza en orden lexicográfico usando la estructura dispersa; CG aplica productos matriz-vector y sus recurrencias estándar. Ninguno reemplaza el algoritmo solicitado por una llamada directa a `pcg`. La división dispersa de MATLAB se usa como referencia algebraica independiente.

La tolerancia del caso base es \(10^{-8}\). Para los ensayos de tiempo, cada configuración tuvo una corrida de calentamiento y tres corridas medidas. Se exportan los tiempos individuales, iteraciones, residuos y estado; las tablas resumen mediana, desviación estándar e IQR. Esta metodología reduce el efecto de una sola medición, pero los tiempos siguen dependiendo del equipo y la carga del sistema.

## Salidas

Las exportaciones se escriben bajo `simulacion/resultados/`: CSV para tablas y tiempos individuales, JSON para manifiestos y metadatos, y MAT para datos MATLAB. El resumen de ejecución está en `simulacion/pruebas/reporte_validacion.json`. El pipeline Python consume estos resultados para producir figuras; MATLAB sigue siendo la fuente de los datos físicos.
