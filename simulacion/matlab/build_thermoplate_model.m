function modelo = build_thermoplate_model(configuracion, nx, ny)
%BUILD_THERMOPLATE_MODEL Compatibilidad con el Live Script archivado.
% Use discretizacion_matriz para el flujo didáctico actual.
    if nargin < 2
        modelo = discretizacion_matriz(configuracion);
    elseif nargin < 3
        modelo = discretizacion_matriz(configuracion, nx);
    else
        modelo = discretizacion_matriz(configuracion, nx, ny);
    end
end
