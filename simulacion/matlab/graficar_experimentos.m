function archivos = graficar_experimentos(ruta_simulacion)
%GRAFICAR_EXPERIMENTOS Crea figuras completas con tablas exportadas por MATLAB.
    set(groot,'defaultFigureColor','white','defaultAxesColor','white', ...
        'defaultAxesXColor','black','defaultAxesYColor','black', ...
        'defaultAxesZColor','black','defaultTextColor','black', ...
        'defaultAxesGridColor',[0.78 0.78 0.78]);
    ruta_resultados=fullfile(ruta_simulacion,'resultados');
    ruta_figuras=fullfile(ruta_resultados,'figuras_matlab');
    if ~isfolder(ruta_figuras), mkdir(ruta_figuras); end
    refinamiento=readtable(fullfile(ruta_resultados,'refinamiento','refinement.csv'));
    sensibilidad=readtable(fullfile(ruta_resultados,'sensibilidad','sensitivity.csv'));
    historial=readtable(fullfile(ruta_resultados,'escenario_base','convergence_history.csv'));
    tiempos=readtable(fullfile(ruta_resultados,'refinamiento','timing_repeats.csv'));

    figura=figure('Name','Refinamiento, tiempos y sensibilidad','Color','w');
    tiledlayout(figura,1,3,'TileSpacing','compact','Padding','compact');
    nexttile; axtoolbar(gca,'Visible','off'); mallas=unique(refinamiento(:,{'nx','ny'}),'rows');
    delta=zeros(height(mallas)-1,1); nodos=zeros(height(mallas)-1,1);
    for j=1:height(mallas)-1
        sel=refinamiento.nx==mallas.nx(j) & refinamiento.ny==mallas.ny(j);
        delta(j)=refinamiento.spatial_max_abs_difference_C(find(sel,1));
        nodos(j)=mallas.nx(j)*mallas.ny(j);
    end
    loglog(nodos,delta,'o-','LineWidth',1.4); grid on
    xlabel('Nodos de la malla'); ylabel('Diferencia maxima (C)'); title('Refinamiento')

    nexttile; axtoolbar(gca,'Visible','off'); hold on
    fuentes=unique(string(sensibilidad.varied_source));
    for k=1:numel(fuentes)
        sel=string(sensibilidad.varied_source)==fuentes(k);
        campo_potencia=sprintf('%s_power_W',char(fuentes(k)));
        valores_potencia=sensibilidad.(campo_potencia);
        plot(valores_potencia(sel),sensibilidad.Tmax_C(sel),'o-', ...
            'LineWidth',1.4,'DisplayName',char(fuentes(k)));
    end
    grid on; xlabel('Potencia (W)'); ylabel('Temperatura maxima (C)');
    title('Sensibilidad'); legend('Location','best')

    nexttile; axtoolbar(gca,'Visible','off'); hold on
    metodos=unique(string(tiempos.solver));
    for k=1:numel(metodos)
        sel=string(tiempos.solver)==metodos(k);
        mallas_t=unique(tiempos.nx(sel)); mediana=zeros(size(mallas_t));
        for j=1:numel(mallas_t)
            mediana(j)=median(tiempos.time_s(sel & tiempos.nx==mallas_t(j)));
        end
        semilogy(mallas_t,mediana,'o-','LineWidth',1.4,'DisplayName',char(metodos(k)));
    end
    grid on; xlabel('Nodos en x'); ylabel('Mediana del tiempo (s)');
    title('Costo de solucion'); legend('Location','best')
    aplicar_estilo_grafico(figura);
    exportgraphics(figura,fullfile(ruta_figuras,'refinamiento_sensibilidad.pdf'), ...
        'ContentType','image','Resolution',220);

    figura_res=figure('Name','Convergencia completa','Color','w'); axtoolbar(gca,'Visible','off'); hold on
    metodos=unique(string(historial.solver));
    for k=1:numel(metodos)
        sel=string(historial.solver)==metodos(k);
        semilogy(historial.iter(sel),historial.rel_residual(sel), ...
            'DisplayName',char(metodos(k)),'LineWidth',1.3);
    end
    grid on; xlabel('Iteracion'); ylabel('Residuo relativo'); title('Convergencia del caso base')
    legend('Location','best')
    aplicar_estilo_grafico(figura_res);
    exportgraphics(figura_res,fullfile(ruta_figuras,'convergencia_completa.pdf'), ...
        'ContentType','image','Resolution',220);
    archivos={fullfile(ruta_figuras,'refinamiento_sensibilidad.pdf'), ...
        fullfile(ruta_figuras,'convergencia_completa.pdf')};
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
