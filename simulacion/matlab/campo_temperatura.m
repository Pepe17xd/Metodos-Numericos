function field = campo_temperatura(model, interior)
%CAMPO_TEMPERATURA Reconstruye el campo nodal completo de temperatura.
    T=model.cfg.Tb*ones(model.ny,model.nx);
    T(2:end-1,2:end-1)=reshape(interior,model.ny-2,model.nx-2);
    q=zeros(model.ny,model.nx);
    q(2:end-1,2:end-1)=model.sourceQ;
    [X,Y]=meshgrid(model.x,model.y);
    field=table(X(:),Y(:),T(:),q(:),'VariableNames', ...
        {'x_m','y_m','T_C','q_W_m3'});
end
