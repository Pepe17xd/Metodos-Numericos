function resultado = ejecutar_proyecto(modo_completo)
%EJECUTAR_PROYECTO Corre el laboratorio rapido y opcionalmente los estudios completos.
% modo_completo=false resuelve el caso base y prepara las figuras esenciales.
% modo_completo=true agrega refinamiento, sensibilidad y pruebas extensas.
    if nargin==0, modo_completo=false; end
    tiempo_total=tic;
    ruta_simulacion=fileparts(mfilename('fullpath'));
    ruta_funcion=fullfile(ruta_simulacion,'matlab');
    addpath(ruta_funcion);

    cfg=planteamiento_modelo();
    modelo=discretizacion_matriz(cfg);
    x0=cfg.initialTemperature*ones(modelo.nUnknowns,1);
    [~,cholFlag]=chol(modelo.A);
    assert(cholFlag==0,'La matriz no paso la comprobacion SPD.');
    assert(modelo.symmetryError<1e-14,'La matriz no es simetrica.');
    assert(abs(modelo.sourcePower.cpu-cfg.cpu.power)<1e-10 && ...
        abs(modelo.sourcePower.gpu-cfg.gpu.power)<1e-10,'La potencia discretizada no se conserva.');
    opciones=struct('tolerance',cfg.tolerance,'maxIterations',5000, ...
        'checkpointFractions',[0.013 0.4 1]);
    fprintf('Malla %d x %d nodos, %d incognitas; dx=%.4g m, dy=%.4g m\n', ...
        modelo.nx,modelo.ny,modelo.nUnknowns,modelo.dx,modelo.dy);
    fprintf('Potencia integrada: CPU %.6g W; GPU %.6g W\n', ...
        modelo.sourcePower.cpu,modelo.sourcePower.gpu);

    % Comparamos la version de referencia con la implementacion eficiente.
    A_demo=[4 -1 -1 0;-1 4 0 -1;-1 0 4 -1;0 -1 -1 4];
    b_demo=[1;2;3;4]; x_demo=zeros(4,1);
    z_prof=gausseidel_profesor(A_demo,b_demo,x_demo,1e-12);
    o_demo=struct('tolerance',1e-12,'maxIterations',1000,'checkpointFractions',[0]);
    r_demo=metodo_gauss_seidel(A_demo,b_demo,x_demo,o_demo);
    fprintf('Equivalencia GS en ejemplo pequeno: %.3g\n', ...
        norm(z_prof(end,1:4)'-r_demo.x,inf));

    % Resolvemos el mismo sistema con ambos metodos y una referencia directa.
    tiempo=tic;
    x_directo=modelo.A\modelo.b;
    gs=metodo_gauss_seidel(modelo.A,modelo.b,x0,opciones);
    cg=metodo_gradiente_conjugado(modelo.A,modelo.b,x0,opciones);
    tiempo_base=toc(tiempo);
    comparacion=comparacion_metodos(modelo,x_directo,gs,cg);
    disp(comparacion)
    assert(strcmp(gs.status,'converged')&&strcmp(cg.status,'converged'), ...
        'Un metodo iterativo no alcanzo la tolerancia. Revise las salidas anteriores.');

    % Dibujamos el campo, la convergencia y las aproximaciones en MATLAB.
    T=crear_matriz_campo(modelo,cg.x);
    figura=figure('Name','Resultados del caso base','Color','w');
    tiledlayout(figura,2,2,'TileSpacing','compact','Padding','compact');
    nexttile; axtoolbar(gca,'Visible','off'); imagesc(modelo.x,modelo.y,T); axis image xy; colorbar
    hold on; dibujar_fuentes(cfg); title('Mapa termico'); xlabel('x (m)'); ylabel('y (m)')
    nexttile; axtoolbar(gca,'Visible','off'); surf(repmat(modelo.x,modelo.ny,1),repmat(modelo.y(:),1,modelo.nx),T,'EdgeColor','none');
    colorbar; title('Superficie estacionaria'); xlabel('x (m)'); ylabel('y (m)'); zlabel('T (C)'); view(40,30)
    nexttile; axtoolbar(gca,'Visible','off'); semilogy(0:gs.iterations,gs.history,'LineWidth',1.3); hold on
    semilogy(0:cg.iterations,cg.history,'LineWidth',1.3); yline(cfg.tolerance,'--');
    grid on; title('Residuo relativo'); xlabel('Iteracion'); ylabel('||b-Ax||/||b||')
    legend('Gauss-Seidel','Gradiente Conjugado','Tolerancia','Location','best')
    nexttile; axtoolbar(gca,'Visible','off'); bar(categorical({'Directa','Gauss-Seidel','Gradiente Conjugado'}), ...
        [max(crear_matriz_campo(modelo,x_directo),[],'all'), ...
         max(crear_matriz_campo(modelo,gs.x),[],'all'),max(T,[],'all')]);
    title('Temperatura maxima (C)'); grid on

    figura_aprox=figure('Name','Aproximaciones iterativas','Color','w');
    tiledlayout(figura_aprox,2,3,'TileSpacing','compact','Padding','compact');
    estados={gs.states,cg.states}; nombres={'Gauss-Seidel','Gradiente Conjugado'};
    for im=1:2
        lista=estados{im}; indice=find([lista.iter]==0,1);
        final=lista(end).iter; medio=find(abs([lista.iter]-final/2)==min(abs([lista.iter]-final/2)),1);
        seleccion=[indice medio numel(lista)];
        for j=1:3
            nexttile; axtoolbar(gca,'Visible','off'); imagesc(modelo.x,modelo.y,crear_matriz_campo(modelo,lista(seleccion(j)).x));
            clim(gca,[cfg.Tb max(T,[],'all')]); axis image xy; colorbar; title(sprintf('%s, k=%d',nombres{im},lista(seleccion(j)).iter));
            xlabel('x (m)'); ylabel('y (m)')
        end
    end
    colormap(figura,turbo); colormap(figura_aprox,turbo)
    % La dependencia lineal permite mostrar una sensibilidad breve en modo rapido.
    figura_sens=figure('Name','Sensibilidad de potencia','Color','w');
    potencias=linspace(0,1.5,7); max_cpu=zeros(size(potencias)); max_gpu=max_cpu;
    for j=1:numel(potencias)
        c=cfg; c.cpu.power=cfg.cpu.power*potencias(j);
        m=discretizacion_matriz(c); max_cpu(j)=max(crear_matriz_campo(m,m.A\m.b),[],'all');
        c=cfg; c.gpu.power=cfg.gpu.power*potencias(j);
        m=discretizacion_matriz(c); max_gpu(j)=max(crear_matriz_campo(m,m.A\m.b),[],'all');
    end
    axtoolbar(gca,'Visible','off'); plot(cfg.cpu.power*potencias,max_cpu,'o-','LineWidth',1.4); hold on
    plot(cfg.gpu.power*potencias,max_gpu,'s-','LineWidth',1.4); grid on
    xlabel('Potencia (W)'); ylabel('Temperatura maxima (C)');
    title('Sensibilidad de la temperatura a la potencia');
    legend('CPU (GPU fija)','GPU (CPU fija)','Location','best')

    carpeta_figuras=fullfile(ruta_simulacion,'resultados','figuras_matlab');
    if ~isfolder(carpeta_figuras), mkdir(carpeta_figuras); end
    aplicar_estilo_grafico(figura);
    aplicar_estilo_grafico(figura_aprox);
    aplicar_estilo_grafico(figura_sens);
    exportgraphics(figura,fullfile(carpeta_figuras,'escenario_base.pdf'),'ContentType','image','Resolution',200);
    exportgraphics(figura_sens,fullfile(carpeta_figuras,'sensibilidad.pdf'),'ContentType','image','Resolution',200);
    exportgraphics(figura_aprox,fullfile(carpeta_figuras,'aproximaciones.pdf'),'ContentType','image','Resolution',200);
    fprintf('Residuo GS %.3g, error de iteracion %.3g (%d iteraciones)\n', ...
        gs.relativeResidual,gs.relativeStepError,gs.iterations);
    fprintf('Residuo CG %.3g (%d iteraciones)\n',cg.relativeResidual,cg.iterations);
    fprintf('Tmax %.6f C; tiempo medido de los solucionadores %.3f s\n',max(T,[],'all'),tiempo_base);

    tiempo_rapido_s=toc(tiempo_total);
    fprintf('Tiempo total del modo rapido, incluidas figuras: %.2f s\n',tiempo_rapido_s);
    resultado=struct('cfg',cfg,'modelo',modelo,'directa',x_directo,'gs',gs,'cg',cg, ...
        'comparacion',comparacion,'tiempoBase_s',tiempo_base, ...
        'tiempoRapido_s',tiempo_rapido_s,'modoCompleto',modo_completo);
    if modo_completo
        fprintf('\nModo completo: se ejecutan refinamiento, sensibilidad extendida y pruebas.\n');
        tiempo_completo=tic;
        resultado.experimentos=experimentos();
        resultado.validacion=validar_simulacion();
        resultado.figurasCompletas=graficar_experimentos(ruta_simulacion);
        resultado.tiempoCompleto_s=toc(tiempo_completo);
        evidencia=struct('experimentsStatus','passed', ...
            'experimentsTime_s',resultado.tiempoCompleto_s, ...
            'experimentScript','simulacion/ejecutar_proyecto.m', ...
            'command','r=ejecutar_proyecto(true)');
        resultado.reporte=generar_reporte_validacion([],[],evidencia);
        disp(resultado.validacion)
    end
end

function dibujar_fuentes(cfg)
% Marca las huellas rectangulares idealizadas de CPU y GPU.
    rectangle('Position',[cfg.cpu.x cfg.cpu.y cfg.cpu.width cfg.cpu.height], ...
        'EdgeColor','w','LineWidth',1.5); text(cfg.cpu.x,cfg.cpu.y,' CPU','Color','w');
    rectangle('Position',[cfg.gpu.x cfg.gpu.y cfg.gpu.width cfg.gpu.height], ...
        'EdgeColor','w','LineWidth',1.5); text(cfg.gpu.x,cfg.gpu.y,' GPU','Color','w');
end


function T=crear_matriz_campo(modelo,x)
% Reconstruye el campo nodal y fija la temperatura del contorno.
    T=modelo.cfg.Tb*ones(modelo.ny,modelo.nx);
    T(2:end-1,2:end-1)=reshape(x,modelo.ny-2,modelo.nx-2);
end








function aplicar_estilo_grafico(figura)
% Usa fondo y texto de alto contraste en cada eje antes de exportar.
    ejes=findall(figura,'Type','axes');
    for j=1:numel(ejes)
        ax=ejes(j); ax.Color='white'; ax.XColor='black'; ax.YColor='black';
        ax.ZColor='black'; ax.GridColor=[0.78 0.78 0.78];
        ax.Title.Color='black'; ax.XLabel.Color='black';
        ax.YLabel.Color='black'; ax.ZLabel.Color='black';
        axtoolbar(ax,'Visible','off');
    end
    barras=findall(figura,'Type','ColorBar');
    for j=1:numel(barras), barras(j).Color='black'; end
end


