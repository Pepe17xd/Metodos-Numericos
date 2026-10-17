%% 1. Presentacion
% Estudiamos la conduccion estacionaria en una placa rectangular idealizada.
% CPU y GPU son regiones con generacion uniforme, no componentes reales.
% La ecuacion es Laplaciano(T) + q'''/k = 0, con T=Tb en el borde.

%% 2. Parametros iniciales
% Dimensiones y propiedades estan en metros, watts y W/(m K).
% La malla base tiene 49 x 33 nodos; la frontera se fija a 25 grados C.
% Las fuentes usan 12 W para CPU y 18 W para GPU.
% Cambie solo modo_completo para ejecutar los estudios largos.
modo_completo=false;

% Buscamos el repositorio desde la carpeta actual y sus carpetas superiores.
ruta_busqueda=pwd;
ruta_laboratorio='';
for nivel=1:8
    if isfile(fullfile(ruta_busqueda,'simulacion','ejecutar_proyecto.m'))
        ruta_laboratorio=fullfile(ruta_busqueda,'simulacion');
        break
    elseif isfile(fullfile(ruta_busqueda,'ejecutar_proyecto.m')) && ...
            isfolder(fullfile(ruta_busqueda,'matlab'))
        ruta_laboratorio=ruta_busqueda;
        break
    end
    ruta_padre=fileparts(ruta_busqueda);
    if strcmp(ruta_padre,ruta_busqueda), break; end
    ruta_busqueda=ruta_padre;
end
assert(~isempty(ruta_laboratorio), ...
    ['No se encontro el repositorio desde la carpeta actual. Abra MATLAB ' ...
     'con la carpeta del repositorio en Current Folder y vuelva a presionar Run.']);
addpath(ruta_laboratorio);
addpath(fullfile(ruta_laboratorio,'matlab'));
cfg=planteamiento_modelo();
disp(cfg)

%% 3. Diferencias Finitas
% Las diferencias centrales aproximan las segundas derivadas.
% En una malla uniforme: T(i,j) es un cuarto de la suma de los cuatro
% vecinos y q'''*h^2/k. Al reunir las ecuaciones interiores se obtiene A*T=b.
% Las temperaturas de borde conocidas se trasladan al vector b.

%% 4. Gauss-Seidel
% Con la convencion del curso A=D-L-U:
% (D-L)*T_nueva = U*T_anterior + b.
% El laboratorio conserva la version literal del profesor en
% matlab/gausseidel_profesor.m. El calculo principal usa la misma recurrencia
% mediante una resolucion triangular dispersa, sin formar una inversa.
% Se registra el error relativo entre iteraciones y el residuo relativo.

%% 5. Gradiente Conjugado
% Para una matriz simetrica definida positiva, CG actualiza x, el residuo r
% y la direccion p con productos dispersos. La implementacion propia esta en
% matlab/metodo_gradiente_conjugado.m; no se usa pcg.

%% 6. Comparacion
% Ambos metodos resuelven la misma matriz, vector y condicion inicial.
% Tambien se compara con el sistema disperso directo de MATLAB.

%% 7. Graficas
% MATLAB crea mapas termicos, superficie, residuos, aproximaciones y
% sensibilidad de temperatura frente a potencia.

%% 8. Refinamiento espacial y validacion
% El modo rapido comprueba simetria, definitud positiva y conservacion de
% potencia. El modo completo ejecuta refinamiento, sensibilidad extendida,
% estadisticas de tiempo y pruebas numericas adicionales.
resultado=ejecutar_proyecto(modo_completo);
