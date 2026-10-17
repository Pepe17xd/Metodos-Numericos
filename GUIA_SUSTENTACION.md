# Guía breve para la sustentación

## ¿Por qué usamos Diferencias Finitas?

La placa es rectangular y la ecuación de conducción es conocida. Las diferencias centrales aproximan las derivadas con temperaturas de nodos vecinos y producen una ecuación sencilla en cada nodo interior.

## ¿De dónde sale la ecuación nodal?

Se sustituyen las aproximaciones centrales de las segundas derivadas en \(\nabla^2T+q'''/k=0\). Si \(\Delta x=\Delta y=h\), al despejar la temperatura se divide entre cuatro toda la suma de los vecinos y la fuente: \(T_{i,j}=\frac14(T_E+T_O+T_N+T_S+q'''h^2/k)\).

## ¿Por qué aparece una matriz?

Cada nodo interior aporta una ecuación que comparte temperaturas con sus vecinos. Al ordenar las incógnitas en un vector, esas ecuaciones forman el sistema disperso \(A T=b\). Las temperaturas conocidas del contorno se trasladan a \(b\).

## ¿Cómo funciona el Gauss–Seidel del profesor?

El laboratorio descompone \(A=D-L-U\) y actualiza usando \((D-L)T^{k+1}=UT^k+b\). La función `gausseidel_profesor.m` conserva literalmente la implementación matricial con inversas para estudiar sistemas pequeños. La función `metodo_gauss_seidel.m` usa la misma recurrencia resolviendo el sistema triangular disperso; no construye las inversas ni la matriz de iteración.

## ¿Qué es una iteración?

Es un paso de actualización de la solución algebraica. No representa un segundo de calentamiento físico. El problema resuelto es estacionario.

## ¿Qué significa el residuo?

El residuo \(r=b-Ax\) mide cuánto incumple una aproximación el sistema discreto. Se usa \(\|r\|_2/\|b\|_2\) para comparar ambos solucionadores con una tolerancia común.

## ¿Cuál es la diferencia entre el error del laboratorio y el residuo?

El error relativo entre aproximaciones de GS, \(\|x^{k+1}-x^k\|_\infty/\|x^{k+1}\|_\infty\), mide cuánto cambió la solución en un paso. El residuo mide el desbalance de \(Ax=b\). Se registran por separado; el residuo común controla la parada.

## ¿Por qué CG requiere una matriz SPD?

Para que \(u^TAv\) defina el producto interno de las direcciones conjugadas y la energía cuadrática tenga un mínimo único. El stencil del laplaciano con temperatura de Dirichlet produce una matriz simétrica definida positiva; el código verifica simetría y factorización de Cholesky.

## ¿Por qué cambia la temperatura con la potencia?

La fuente térmica aumenta con la potencia. Como el sistema es lineal con propiedades y fronteras fijas, las temperaturas cambian con la generación. El máximo global puede moverse entre CPU y GPU y no tiene por qué variar con una sola pendiente.

## ¿Por qué refinamos la malla?

Al aumentar el número de nodos observamos cuánto cambia la solución espacial. La diferencia entre mallas sucesivas es un indicador de refinamiento, no un error exacto porque no se cuenta con una solución analítica del modelo completo.

## ¿Por qué el modelo no describe exactamente una laptop?

Representa una placa homogénea 2D, estacionaria y con conductividad constante, fuentes rectangulares uniformes y todo el borde fijado a una temperatura. No modela ventiladores, disipadores, convección, radiación, contactos, transitorios ni capas reales; las temperaturas son resultados del modelo idealizado.

## Resultados de referencia recalculados

En MATLAB R2025a Update 1, para la malla \(49\times33\), el máximo directo fue \(77.759239\ ^\circ\mathrm C\). GS convergió en 2145 iteraciones y CG en 127, ambos con residuo relativo menor que \(10^{-8}\). La potencia discreta fue 12 W para CPU y 18 W para GPU. Son resultados simulados, no mediciones experimentales.
