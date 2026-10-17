function [z] = gausseidel_profesor(A,b,x0,Tol)
%GAUSSEIDEL_PROFESOR Implementacion literal entregada en el laboratorio.
    D = diag(diag(A));
    L = -tril(A,-1);
    U = -triu(A,1);
    Tgs = inv(D-L)*U;
    cgs = inv(D-L)*b;
    error = 1;
    z = [x0' error];
    while error > Tol
        x1 = Tgs*x0 + cgs;
        error = norm(x1-x0,'inf')/norm(x1,'inf');
        z = [z; x1' error];
        x0 = x1;
    end
end
