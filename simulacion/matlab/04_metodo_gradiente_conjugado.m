function out = metodo_gradiente_conjugado(A,b,x0,opts)
%METODO_GRADIENTE_CONJUGADO Gradiente Conjugado explícito sin precondicionar.
%   Solves the sparse symmetric positive-definite system without calling
%   pcg. opts may contain tolerance, maxIterations, and checkpointFractions.
%   The reported residual is ||b-A*x||_2/||b||_2, using denominator one
%   only for the zero right-hand side.
%   For f(x)=0.5*x'*A*x-b'*x, r=b-A*x is minus the gradient. The search
%   directions combine the current residual with the preceding direction.
% The routine checks positive curvature and reports breakdown rather than
% silently continuing when assumptions or floating-point arithmetic fail.
    if nargin<4, opts=struct(); end
    opts=solver_options(opts); n=numel(b); x=x0(:);
    assert(isequal(size(A),[n n]) && numel(x)==n,'A, b, and x0 sizes do not match.');
    bn=norm(b); if bn==0,bn=1;end
    r=b-A*x; p=r; rr=r'*r; rel=sqrt(rr)/bn;
    history=nan(opts.maxIterations+1,1); history(1)=rel;
    checkpoints=iteration_checkpoints(opts);
    states=struct('iter',0,'x',x); status='max_iterations';
    if rel<=opts.tolerance
        out=struct('solver','Conjugate Gradient','x',x,'iterations',0, ...
            'relativeResidual',rel,'time_s',0,'status','converged', ...
            'history',history(1),'states',states,'tolerance',opts.tolerance, ...
            'relativeStepError',NaN); return;
    end
    t=tic;
    for iter=1:opts.maxIterations
        Ap=A*p; curvature=p'*Ap;
        if ~isfinite(curvature) || curvature<=0, status='nonpositive_curvature'; break; end
        alpha=rr/curvature; x=x+alpha*p; r=r-alpha*Ap;
        rrNew=r'*r; rel=sqrt(max(rrNew,0))/bn; history(iter+1)=rel;
        if any(checkpoints==iter), states(end+1)=struct('iter',iter,'x',x); end %#ok<AGROW>
        if ~isfinite(rel) || any(~isfinite(x)), status='nonfinite'; break; end
        if rel<=opts.tolerance, status='converged'; rr=rrNew; break; end
        if rr==0, status='zero_residual_breakdown'; break; end
        beta=rrNew/rr; p=r+beta*p; rr=rrNew;
    end
    elapsed=toc(t); nDone=iter; history=history(1:nDone+1);
    relExact=norm(b-A*x)/bn;
    if states(end).iter~=nDone, states(end+1)=struct('iter',nDone,'x',x); end
    out=struct('solver','Conjugate Gradient','x',x,'iterations',nDone, ...
        'relativeResidual',relExact,'time_s',elapsed,'status',status, ...
        'history',history,'states',states,'tolerance',opts.tolerance, ...
        'relativeStepError',NaN);
end

function o=solver_options(o)
    defaults=struct('tolerance',1e-8,'maxIterations',200000, ...
        'checkpointFractions',[0.01 0.1 0.5 1]);
    f=fieldnames(defaults); for j=1:numel(f), if ~isfield(o,f{j}),o.(f{j})=defaults.(f{j});end,end
    validateattributes(o.tolerance,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(o.maxIterations,{'numeric'},{'scalar','real','finite','integer','>=',1});
    validateattributes(o.checkpointFractions,{'numeric'},{'vector','real','finite','>=',0,'<=',1});
end

function checkpoints=iteration_checkpoints(opts)
% Include zero and the final iteration so each run has endpoint states.
    checkpoints=unique([0 round(opts.checkpointFractions(:)'*opts.maxIterations) ...
        opts.maxIterations]);
end
