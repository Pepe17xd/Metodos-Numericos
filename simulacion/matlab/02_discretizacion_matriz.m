function model = discretizacion_matriz(cfg, nx, ny)
%DISCRETIZACION_MATRIZ Ensambla el sistema disperso de diferencias finitas.
%   MODEL = BUILD_THERMOPLATE_MODEL(CFG,NX,NY) returns A and b for
%   -laplacian(T)=q'''/k with fixed boundary temperature CFG.Tb. Lengths
%   are metres, conductivity is W/(m K), source power is watts, and
%   temperature is reported in degC. NX and NY count nodes including edges.
% Unknowns are interior nodes in MATLAB column-major order (y varies fastest).
% Each interior-node source is averaged over its dual control volume. The
% rectangle/control-volume overlap fractions preserve integrated power even
% when configurable source edges fall between grid nodes.
    if nargin < 2 || isempty(nx), nx = cfg.baselineGrid(1); end
    if nargin < 3 || isempty(ny), ny = cfg.baselineGrid(2); end
    validateattributes(nx,{'numeric'},{'scalar','integer','>=',3});
    validateattributes(ny,{'numeric'},{'scalar','integer','>=',3});
    validateattributes(cfg.Lx,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(cfg.Ly,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(cfg.k,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(cfg.thickness,{'numeric'},{'scalar','real','finite','positive'});
    x = linspace(0,cfg.Lx,nx); y = linspace(0,cfg.Ly,ny);
    dx=x(2)-x(1); dy=y(2)-y(1); nxi=nx-2; nyi=ny-2;
    xlo=x(2:end-1)-dx/2; xhi=x(2:end-1)+dx/2;
    ylo=y(2:end-1)-dy/2; yhi=y(2:end-1)+dy/2;
    q = zeros(nyi,nxi); sourcePower = struct();
    names = {'cpu','gpu'};
    for s=1:numel(names)
        name=names{s}; r=cfg.(name);
        assert(r.x>=xlo(1) && r.x+r.width<=xhi(end) && ...
            r.y>=ylo(1) && r.y+r.height<=yhi(end), ...
            'Source rectangles must lie within active interior-node control volumes.');
        qvol = r.power/(r.width*r.height*cfg.thickness);
        overlapX=max(0,min(xhi,r.x+r.width)-max(xlo,r.x));
        overlapY=max(0,min(yhi,r.y+r.height)-max(ylo,r.y));
        overlapArea=overlapY(:)*overlapX(:)';
        quadratureFactor=overlapArea/(dx*dy);
        q=q+qvol*quadratureFactor;
        sourcePower.(name)=sum(qvol*overlapArea,'all')*cfg.thickness;
    end
    % A*T=b where A approximates -Laplacian. Dirichlet contributions are
    % moved to b; the RHS source term is q'''/k.
    % Kronecker assembly follows MATLAB's column-major vectorization:
    % reshape(T,nyi,nxi) makes y the fast index in the unknown vector.
    Tx = spdiags([-ones(nxi,1)/dx^2, 2*ones(nxi,1)/dx^2, ...
        -ones(nxi,1)/dx^2],[-1 0 1],nxi,nxi);
    Ty = spdiags([-ones(nyi,1)/dy^2, 2*ones(nyi,1)/dy^2, ...
        -ones(nyi,1)/dy^2],[-1 0 1],nyi,nyi);
    A = kron(speye(nxi),Ty) + kron(Tx,speye(nyi));
    rhs = q(:)/cfg.k;
    boundary = zeros(nyi,nxi);
    boundary(:,1)=boundary(:,1)+cfg.Tb/dx^2;
    boundary(:,end)=boundary(:,end)+cfg.Tb/dx^2;
    boundary(1,:)=boundary(1,:)+cfg.Tb/dy^2;
    boundary(end,:)=boundary(end,:)+cfg.Tb/dy^2;
    b=rhs+boundary(:);
    model=struct('cfg',cfg,'nx',nx,'ny',ny,'x',x,'y',y,'dx',dx,'dy',dy, ...
        'nUnknowns',nxi*nyi,'A',A,'b',b,'sourceQ',q,'sourceTerm',rhs, ...
        'sourcePower',sourcePower,'symmetryError',norm(A-A','fro')/max(1,norm(A,'fro')));
end
