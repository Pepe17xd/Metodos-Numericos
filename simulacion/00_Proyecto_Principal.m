%% Modelación térmica de una placa electrónica idealizada
% Curso CC2104, Métodos Numéricos. Las regiones CPU y GPU son fuentes
% rectangulares idealizadas. Sus potencias son parámetros de simulación,
% no mediciones de una computadora.

clear; clc;
% Las rutas se obtienen desde este archivo: no dependen de la carpeta actual.
ruta_script = mfilename('fullpath');
if isempty(ruta_script)
    ruta_simulacion = fullfile(pwd,'simulacion');
else
    ruta_simulacion = fileparts(ruta_script);
end
ruta_matlab = fullfile(ruta_simulacion,'matlab');
assert(isfile(fullfile(ruta_matlab,'planteamiento_modelo.m')), ...
    'Abra 00_Proyecto_Principal.mlx desde una copia completa del repositorio.');
addpath(ruta_matlab);

%% 1. Problema físico y parámetros
% Se estudia una placa homogénea e isótropa en estado estacionario. El calor
% se propaga por conducción en 2-D. Se fija T=Tb en todo el contorno y se
% omiten convección, radiación, contacto térmico, capas y transitorios.
%
% Ecuación: ∂²T/∂x² + ∂²T/∂y² + q'''/k = 0 en Ω; T=Tb en ∂Ω.
%
% CPU y GPU se representan como dos regiones con generación volumétrica
% uniforme. modelo_base.json contiene dimensiones, conductividad, espesor,
% temperatura de borde, huellas, potencias y mallas del estudio.
cfg = planteamiento_modelo();
disp(cfg)

%% 2. Diferencias finitas y ensamblaje disperso
% Para espaciamiento uniforme h, la diferencia central de segundo orden da
% T(i,j)=1/4*[T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)+q'''h²/k].
% El factor 1/4 multiplica toda la suma. Para Δx y Δy distintos se usan los
% pesos 1/Δx² y 1/Δy². El ensamblador crea el stencil de cinco puntos para
% -∇²T=q'''/k. Las incógnitas son los nodos interiores; Dirichlet pasa a b.
model = discretizacion_matriz(cfg);
fprintf('Malla: %d x %d nodos; incógnitas: %d; dx=%.4g m; dy=%.4g m\n', ...
    model.nx,model.ny,model.nUnknowns,model.dx,model.dy);
fprintf('Potencia discretizada: CPU %.12g W; GPU %.12g W\n', ...
    model.sourcePower.cpu,model.sourcePower.gpu);

%% 3. Ejemplo explícito reducido de A*T=b
% Ejemplo de 4 incógnitas en una malla 4 x 4 con Δx=Δy=1 y Tb=25.
% El orden vectorial es column-major, con y variando más rápido.
% El stencil produce la matriz y las contribuciones de borde suman 50 K
% por fila. La parte restante de b representa q'''/k en este ejemplo.
A_demo = [ 4 -1 -1  0; ...
          -1  4  0 -1; ...
          -1  0  4 -1; ...
           0 -1 -1  4];
source_demo = [0;2;1;3];
b_demo = 50*ones(4,1) + source_demo;
disp('Ejemplo reducido: A_demo ='); disp(A_demo)
disp('Ejemplo reducido: b_demo ='); disp(b_demo)

%% 4. Algoritmo de Gauss–Seidel
% Se usa la convención del laboratorio: A=D-L-U,
% D=diag(diag(A)), L=-tril(A,-1), U=-triu(A,1).
% (D-L)T^(k+1)=U*T^k+b; T^(k+1)=T_GS*T^k+c_GS,
% T_GS=(D-L)^(-1)U y c_GS=(D-L)^(-1)b.
% La inversa es notación matemática: en el sistema disperso se resuelve un
% sistema triangular. Componente a componente se reutilizan los valores
% nuevos para j<i y los anteriores para j>i.
%
% La convergencia de una iteración estacionaria depende de ρ(T_GS)<1.
% La dominancia diagonal estricta es suficiente, pero no necesaria; no
% inferimos divergencia solo por su ausencia. Aquí A es simétrica definida
% positiva por la discretización de Dirichlet y se comprueba con chol.

%% 5. Gradiente Conjugado: energía y recurrencias
% CG requiere A simétrica definida positiva. Entonces
% f(T)=1/2*T'*A*T-b'*T es estrictamente convexa y su minimizador único
% satisface A*T=b. ∇f=A*T-b y r=b-A*T=-∇f es el residuo.
% El descenso por gradiente usa solo el gradiente actual; CG combina el
% residuo con la dirección previa para producir direcciones A-conjugadas.
% Se implementan explícitamente las recurrencias:
% alpha_k=(r_k'*r_k)/(p_k'*A*p_k), T_(k+1)=T_k+alpha_k*p_k,
% r_(k+1)=r_k-alpha_k*A*p_k,
% beta_k=(r_(k+1)'*r_(k+1))/(r_k'*r_k),
% p_(k+1)=r_(k+1)+beta_k*p_k.
% No se llama a pcg. La combinación de direcciones evita repetir, en la
% métrica A, correcciones ya realizadas; en los casos medidos esto redujo
% el conteo frente a GS, sin implicar una ventaja universal.

%% 6. Criterios de convergencia
% El criterio de parada común es el residuo relativo euclídeo
% rho_k=norm(b-A*T_k,2)/norm(b,2), con tolerancia configurada 1e-8.
% Para GS se registra además E_k=norm(T_(k+1)-T_k,inf)/norm(T_(k+1),inf).
% E_k mide el cambio entre aproximaciones, no directamente el desbalance
% A*T=b; se informa como complemento y no modifica la parada histórica.
x0 = cfg.initialTemperature*ones(model.nUnknowns,1);
opts = struct('tolerance',cfg.tolerance,'maxIterations',cfg.maxIterations, ...
    'checkpointFractions',cfg.checkpointFractions);
[~,cholFlag] = chol(model.A);
assert(cholFlag==0,'La matriz ensamblada no pasó la comprobación SPD.');

%% 7. Solución y comparación en el caso base
% La división dispersa aporta una referencia independiente. GS y CG reciben
% la misma matriz, vector y condición inicial.
baseTimer = tic;
x_direct = model.A\model.b;
gs = metodo_gauss_seidel(model.A,model.b,x0,opts);
cg = metodo_gradiente_conjugado(model.A,model.b,x0,opts);
comparison = comparacion_metodos(model,x_direct,gs,cg);
disp(comparison)
baseElapsed_s = toc(baseTimer);

%% 8. Exportación de resultados para las figuras Python
% El orquestador ejecuta el escenario base, refinamiento, sensibilidad y
% pruebas de tolerancia. Es MATLAB quien resuelve y exporta CSV/JSON/MAT;
% visualizar.py solo lee esos datos para crear los gráficos.
fprintf(['Las exportaciones aparecerán en simulacion/resultados/. ' ...
    'El flujo completo se inicia en la sección siguiente.\n']);
experimentTimer = tic;
experimentReport = experimentos();
experimentsElapsed_s = toc(experimentTimer);
executionEvidence = struct('baseStatus','passed','baseTime_s',baseElapsed_s, ...
    'experimentsStatus','passed','experimentsTime_s',experimentsElapsed_s, ...
    'experimentScript','simulacion/00_Proyecto_Principal.mlx', ...
    'command','matlab.internal.liveeditor.executeAndSave(''simulacion/00_Proyecto_Principal.mlx'')');
validationReport = generar_reporte_validacion([],[],executionEvidence);
fprintf('Exportaciones completas: %s\n',experimentReport.baseDirectory);
fprintf('Validación: %s; tiempo del flujo completo: %.2f s\n', ...
    validationReport.validation_status,experimentsElapsed_s);
