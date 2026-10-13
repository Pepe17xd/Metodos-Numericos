function reporte = run_thermoplate_experiments(ruta_configuracion, ruta_salida, opciones_tiempo)
%RUN_THERMOPLATE_EXPERIMENTS Compatibilidad; use ejecutar_experimentos.
    if nargin < 1
        reporte = ejecutar_experimentos();
    elseif nargin < 2
        reporte = ejecutar_experimentos(ruta_configuracion);
    elseif nargin < 3
        reporte = ejecutar_experimentos(ruta_configuracion, ruta_salida);
    else
        reporte = ejecutar_experimentos(ruta_configuracion, ruta_salida, opciones_tiempo);
    end
end
