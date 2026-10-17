# Proyecto de Métodos Numéricos: placa térmica idealizada

El proyecto calcula la temperatura estacionaria en una placa rectangular homogénea con dos regiones de generación de calor asociadas con CPU y GPU. Es un modelo matemático didáctico; sus parámetros no son mediciones de una laptop.

## Métodos numéricos

- **Diferencias Finitas:** aproxima la ecuación de conducción y construye el sistema disperso `A*T = b`.
- **Gauss–Seidel:** actualiza las temperaturas en orden. Se conserva la versión literal vista en clase y se usa una versión eficiente con la misma recurrencia.
- **Gradiente Conjugado:** implementa explícitamente las iteraciones para la matriz simétrica definida positiva. No usa `pcg`.

La solución dispersa directa de MATLAB sirve como referencia numérica. El error relativo entre iteraciones de Gauss–Seidel y el residuo relativo común se reportan por separado.

## Ejecutar

1. Abra `simulacion/00_Proyecto_Principal.mlx` en MATLAB.
2. Presione **Run**. El Live Script ejecuta el modo rápido, resuelve el escenario base y muestra los resultados y gráficas.

El Live Script localiza el proyecto si la carpeta actual de MATLAB está en la raíz del repositorio o dentro de sus carpetas. Si MATLAB se abrió desde otra ubicación, abra primero la carpeta del repositorio. No necesita instalar paquetes ni generar archivos previos.

Para ejecutar sin Live Editor, abra `simulacion/ejecutar_proyecto.m` y presione **Run**, o invoque `ejecutar_proyecto(false)` desde MATLAB cuando la carpeta `simulacion` esté en el path.

## Modo completo

Al inicio de `00_Proyecto_Principal.mlx`, cambie una sola vez:

```matlab
modo_completo = true;
```

El modo completo repite las cuatro mallas, el análisis de sensibilidad, las repeticiones para tiempos y las pruebas numéricas. Regenera tablas, metadatos y figuras dentro de `simulacion/resultados/`. Los tiempos se miden en cada ejecución; los CSV históricos no se presentan como tiempos nuevos.

## Resultados

El modo rápido muestra la malla, potencias integradas, temperatura máxima, comparación de solucionadores, error de iteraciones de Gauss–Seidel, residuos, aproximaciones inicial/intermedia/final y figuras. Los PDF creados por MATLAB quedan en `simulacion/resultados/figuras_matlab/`.

El modo completo añade el refinamiento espacial, la sensibilidad extendida, pruebas de convergencia y tiempos medidos. Las exportaciones históricas existentes se conservan para comparación; no sustituyen una ejecución nueva.

## Compatibilidad comprobada

El flujo se ejecutó en **MATLAB R2025a Update 1 (25.1.0.2973910)** en Windows. No se ha comprobado en otras versiones o sistemas operativos.

## Archivos necesarios

- `simulacion/00_Proyecto_Principal.mlx`: laboratorio principal.
- `simulacion/ejecutar_proyecto.m`: cálculo compartido y modo sin Live Editor.
- `simulacion/matlab/`: discretización, solucionadores y verificaciones.
- `simulacion/resultados/`: las ejecuciones nuevas se calculan sin depender de los resultados guardados; las salidas se conservan como evidencia histórica.

No se requiere Python, JSON de configuración, toolbox adicional ni ruta absoluta.

## Limitaciones

La placa es 2D, homogénea e isotrópica, estacionaria y de conductividad constante. Se prescribe la temperatura en todo el borde y se omiten convección, radiación, ventiladores, contacto, transitorios y capas de una placa real. Por ello, las temperaturas calculadas no predicen directamente el funcionamiento de una computadora.


