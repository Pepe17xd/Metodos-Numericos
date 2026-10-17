# Funciones MATLAB

Las funciones de esta carpeta implementan la malla, la matriz dispersa, los solucionadores y las verificaciones del proyecto. Los parámetros están definidos en `configuracion_modelo.m`; no hace falta un archivo de configuración externo.

- `discretizacion_matriz.m`: crea la malla, integra CPU y GPU y ensambla `A*T=b` con fronteras de Dirichlet.
- `gausseidel_profesor.m`: copia literal de la función mostrada en el laboratorio del curso; se usa solo como referencia pequeña.
- `metodo_gauss_seidel.m`: misma recurrencia con resolución triangular dispersa, residuo común, error de iteración y límite de iteraciones.
- `metodo_gradiente_conjugado.m`: recurrencias propias de CG, con revisión de curvatura y residuo.
- `experimentos.m`: modo completo, con refinamiento, sensibilidad, repeticiones y metadatos.
- `validar_simulacion.m`: comprueba fuentes, frontera, matriz SPD, convergencia y referencia directa.

`ejecutar_proyecto(false)` realiza el modo rápido. `ejecutar_proyecto(true)` ejecuta también el modo completo. Las figuras se generan desde MATLAB.
