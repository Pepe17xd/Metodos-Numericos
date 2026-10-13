function m = metricas_temperatura(model,x)
%METRICAS_TEMPERATURA Calcula máximo térmico y residuo algebraico.
    T=model.cfg.Tb*ones(model.ny,model.nx);
    T(2:end-1,2:end-1)=reshape(x,model.ny-2,model.nx-2);
    maxT=max(T,[],'all'); [iy,ix]=find(T==maxT);
    m=struct('Tmax_C',maxT,'hotspot_x_m',model.x(ix(1)), ...
        'hotspot_y_m',model.y(iy(1)), ...
        'relative_residual',norm(model.b-model.A*x)/max(norm(model.b),1), ...
        'finite',all(isfinite(T),'all'));
end
