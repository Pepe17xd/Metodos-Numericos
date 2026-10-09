# Visualización científica THERMOPLATE

`visualizar.py` genera figuras estáticas a partir de archivos exportados por
MATLAB. Python no ensambla matrices, resuelve el sistema térmico ni modifica los
resultados. Las figuras se guardan en PNG a 300 dpi y PDF vectorial cuando la
entrada correspondiente está disponible.

## Instalación y ejecución

Requiere Python 3.10 o posterior, NumPy, Pandas y Matplotlib. Desde la raíz
del proyecto, instale dependencias y genere las figuras a partir de las
exportaciones MATLAB existentes:

```powershell
python -m pip install -r simulacion/python/requirements.txt
python simulacion/python/visualizar.py --strict --results simulacion/resultados --output simulacion/figuras
```

También puede ejecutar el script desde esta carpeta:

```powershell
python -m pip install -r requirements.txt
python visualizar.py --strict --results ..\resultados --output ..\figuras
```

Sin argumentos, el lector busca `simulacion/resultados/` y escribe en
`simulacion/python/figuras/`.

Las pruebas ligeras del lector usan tablas temporales y no generan figuras ni
resultados científicos:

```powershell
python -B -m unittest discover -s simulacion/python/tests -v
```

`--strict` exige los cinco directorios de escenario y termina con error si
faltan entradas o metadatos declarados al cargar los manifiestos. En el modo normal se crean solo las
figuras sustentadas por datos existentes y se informa qué quedó omitido. Si no
hay entradas utilizables, el programa termina con código 2 y no crea figuras
de muestra.

## Intercambio esperado

MATLAB exporta bajo `simulacion/resultados/` estos directorios:

```text
escenario_base/
convergencia/
refinamiento/
sensibilidad/
estres_numerico/
```

Cada directorio puede incluir un `manifest.json` con las rutas relativas de
sus tablas. El programa también busca archivos CSV, JSON y MAT de forma
recursiva cuando no hay manifiesto, usando el nombre de archivo como nombre
lógico del dataset. Ejemplo mínimo:

```json
{
  "experiment": "escenario_base",
  "datasets": {
    "field": "campo_termico.csv"
  },
  "parameters": {
    "geometry": {
      "Lx_m": 0.24,
      "Ly_m": 0.16,
      "sources": [
        {"name": "CPU", "x_m": 0.04, "y_m": 0.08, "width_m": 0.04, "height_m": 0.04},
        {"name": "GPU", "x_m": 0.16, "y_m": 0.08, "width_m": 0.04, "height_m": 0.04}
      ]
    }
  }
}
```

Los nombres de columnas distinguen mayúsculas y minúsculas solo en los datos
de salida, no durante su detección. El esquema canónico es:

| Escenario / dataset | Columnas requeridas |
|---|---|
| `escenario_base/field` | `x_m,y_m,T_C` |
| `convergencia/history` | `iter,solver,rel_residual` |
| `convergencia/states` | `iter,solver,x_m,y_m,T_C` para evolución espacial |
| `refinamiento/mesh` | `nx,ny,unknowns,solver,time_s,rel_residual` |
| `sensibilidad/sensitivity` | `cpu_power_W,gpu_power_W,Tmax_C` |
| `estres_numerico/stress` | `case,tol,solver,status,iterations,residual,time_s` |

El manifiesto debe declarar los datasets bajo `datasets`, por ejemplo
`"history": "historial.csv"`. Para el esquema geométrico, incluya
`parameters.geometry` o `metadata.geometry` con `Lx_m`, `Ly_m` y una lista
`sources`; cada fuente tiene `name`, `x_m`, `y_m`, `width_m` y `height_m`.
La tabla `field` puede exportarse como filas largas (`x_m,y_m,T_C`) o como una
grilla JSON en `field_grid` con vectores `x`, `y` y matriz `T` de forma
`(len(y),len(x))`.

Las tablas JSON pueden ser listas de filas, un objeto con `data`/`rows`, o
columnas representadas como listas. CSV es el formato recomendado. MAT se
admite para tablas simples cuando SciPy está instalado; se recomienda exportar
CSV/JSON desde MATLAB para evitar limitaciones de estructuras MAT.

## Figuras

1. Esquema geométrico y fuentes (requiere geometría en metadatos).
2. Mapa térmico 2D y 3. superficie 3D (requieren `field`).
4. Evolución de estados iterativos (requiere `states`).
5. Comparación visual Gauss–Seidel / Gradiente Conjugado y 6. curvas de
   residuo relativo (requieren `history`).
7. Costo computacional por refinamiento (`mesh`), 8. sensibilidad térmica
   (`sensitivity`) y 9. resultados de estrés numérico (`stress`).

Si una exportación no contiene una tabla, columnas o metadatos necesarios, se
omite solo la figura afectada y se muestra un aviso. Iteraciones representan
progreso algebraico; no son tiempo físico de calentamiento.
