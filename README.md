# Modelación numérica de la distribución térmica en una placa idealizada

**Universidad de Ingeniería y Tecnología (UTEC)**  
**Curso:** CC2104 – Métodos Numéricos  
**Herramienta de cálculo y visualización:** MATLAB

Este proyecto estudia la temperatura estacionaria de una placa rectangular homogénea con dos regiones de generación de calor que representan CPU y GPU. Es una idealización matemática; los parámetros no son mediciones de una laptop.

## Métodos numéricos

- **Diferencias Finitas:** aproxima la ecuación de conducción y construye el sistema disperso `A*T = b`.
- **Gauss–Seidel:** conserva la implementación literal del curso como referencia y resuelve mallas grandes con una variante eficiente basada en la misma recurrencia.
- **Gradiente Conjugado:** implementa explícitamente las iteraciones para la matriz simétrica definida positiva; no usa `pcg`.

La solución dispersa directa de MATLAB se usa como referencia numérica. El error relativo entre iteraciones de Gauss–Seidel y el residuo relativo común se informan por separado. MATLAB genera directamente las gráficas. El proyecto ya no requiere Python, sus paquetes ni scripts auxiliares.

## Inicio rápido: descargar y ejecutar

No se necesita Git para usar el proyecto:

1. Abre el [repositorio en GitHub](https://github.com/Pepe17xd/Metodos-Numericos) y selecciona **Code → Download ZIP**.
2. Extrae el ZIP y abre la carpeta del proyecto en MATLAB.
3. Abre `simulacion/00_Proyecto_Principal.mlx` y pulsa **Run**.

El modo rápido resuelve el escenario base y muestra las métricas y gráficas esenciales. No necesita archivos de resultados previos, configuración JSON, instalación de paquetes ni toolboxes adicionales. Al ejecutar el Live Script desde el Editor, este localiza sus funciones junto al archivo aunque la carpeta actual de MATLAB sea otra.

Para ejecutar sin Live Editor, abre `simulacion/ejecutar_proyecto.m` y pulsa **Run**.

## Modo completo

Al inicio de `00_Proyecto_Principal.mlx`, cambia una sola variable:

```matlab
modo_completo = true;
```

Este modo repite las cuatro mallas, el análisis de sensibilidad, las repeticiones de tiempos y las pruebas numéricas. Regenera resultados, metadatos, tablas y figuras en `simulacion/resultados/`. Los tiempos se vuelven a medir; los CSV guardados no sustituyen una ejecución nueva.

## Resultados verificados del escenario base

Para la malla de 49 × 33 nodos (1457 incógnitas interiores), con 12 W en CPU y 18 W en GPU:

| Indicador | Resultado |
|---|---:|
| Temperatura máxima de referencia directa | 77.759239 °C |
| Gauss–Seidel | 2145 iteraciones; residuo relativo `9.9420e-9` |
| Gradiente Conjugado | 127 iteraciones; residuo relativo `9.2837e-9` |
| Tolerancia de residuo | `1e-8` |

Estos resultados se recalcularon en MATLAB R2025a Update 1 (25.1.0.2973910) en Windows. No se ha comprobado la ejecución en otras versiones o sistemas operativos. Los tiempos dependen del equipo.

## Organización del repositorio

- `simulacion/00_Proyecto_Principal.mlx`: laboratorio principal, ejecutable con un Run.
- `simulacion/ejecutar_proyecto.m`: cálculo compartido por el Live Script y entrada sin Live Editor.
- `simulacion/matlab/`: discretización, Gauss–Seidel, Gradiente Conjugado, experimentos y validaciones.
- `simulacion/resultados/`: evidencia numérica, datos y figuras generadas desde MATLAB.
- `GUIA_SUSTENTACION.md`: preguntas frecuentes para explicar el modelo y los métodos.

## Alcance y limitaciones

El modelo es bidimensional, estacionario, homogéneo e isotrópico, con conductividad constante y temperatura prescrita en todo el borde. No incorpora convección, radiación, ventiladores, disipadores, contactos, transitorios ni capas de una placa real. Por ello, sus temperaturas describen el modelo matemático y no predicen directamente las de una computadora.

## Artículo científico de referencia

Gonçalves de Brito dos Santos, V., & Gomes dos Anjos, P. T. (2022). Finite difference method applied in two-dimensional heat conduction problem in the permanent regime in rectangular coordinates. *Advances in Pure Mathematics, 12*, 505–518. https://doi.org/10.4236/apm.2022.129038

El artículo aporta el antecedente de conducción y diferencias finitas; la placa idealizada con fuentes CPU/GPU y la comparación de solucionadores son la adaptación del grupo.
