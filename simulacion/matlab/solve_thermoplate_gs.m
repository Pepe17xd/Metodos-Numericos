function out = solve_thermoplate_gs(A,b,x0,opts)
%SOLVE_THERMOPLATE_GS Explicit lexicographic Gauss--Seidel iteration.
%   Updates x from (D+L)x_new = b-U*x_old, where A=D+L+U. For the
%   five-point plate matrix the lower-triangular solve is sparse. opts may
%   contain tolerance, maxIterations, residualCheckInterval, and
%   checkpointFractions. The returned residual is ||b-A*x||_2/||b||_2
%   (with denominator one only when b is exactly zero).
    if nargin<4, opts=struct(); end
    opts=solver_options(opts); n=numel(b); x=x0(:);
    assert(isequal(size(A),[n n]) && numel(x)==n,'A, b, and x0 sizes do not match.');
    d=diag(A); assert(all(isfinite(d))&&all(d~=0),'A has a zero/nonfinite diagonal.');
    M=spdiags(d,0,n,n)+tril(A,-1); U=triu(A,1);
    bn=norm(b); if bn==0, bn=1; end
    checkpoints=iteration_checkpoints(opts);
    history=nan(opts.maxIterations+1,1); history(1)=norm(b-A*x)/bn;
    states=struct('iter',{},'x',{}); states(1)=struct('iter',0,'x',x);
    if history(1)<=opts.tolerance, out=finish('converged',0,0); return; end
    status='max_iterations'; t=tic;
    for iter=1:opts.maxIterations
        x=M\(b-U*x);
        if mod(iter,opts.residualCheckInterval)==0 || iter==opts.maxIterations || any(checkpoints==iter)
            rel=norm(b-A*x)/bn;
            history(iter+1)=rel;
            if any(checkpoints==iter), states(end+1)=struct('iter',iter,'x',x); end %#ok<AGROW>
            if ~isfinite(rel) || any(~isfinite(x)), status='nonfinite'; break; end
            if rel<=opts.tolerance, status='converged'; break; end
        end
    end
    elapsed=toc(t); nDone=iter;
    history=history(1:nDone+1);
    if states(end).iter~=nDone, states(end+1)=struct('iter',nDone,'x',x); end
    out=struct('solver','Gauss-Seidel','x',x,'iterations',nDone, ...
        'relativeResidual',norm(b-A*x)/bn,'time_s',elapsed,'status',status, ...
        'history',history,'states',states,'tolerance',opts.tolerance);
    function o=finish(st,it,tm)
        o=struct('solver','Gauss-Seidel','x',x,'iterations',it, ...
            'relativeResidual',history(1),'time_s',tm,'status',st, ...
            'history',history(1),'states',states,'tolerance',opts.tolerance);
    end
end

function o=solver_options(o)
    defaults=struct('tolerance',1e-8,'maxIterations',200000, ...
        'checkpointFractions',[0.01 0.1 0.5 1],'residualCheckInterval',1);
    f=fieldnames(defaults); for j=1:numel(f), if ~isfield(o,f{j}),o.(f{j})=defaults.(f{j});end,end
    validateattributes(o.tolerance,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(o.maxIterations,{'numeric'},{'scalar','real','finite','integer','>=',1});
    validateattributes(o.residualCheckInterval,{'numeric'},{'scalar','real','finite','integer','>=',1});
    validateattributes(o.checkpointFractions,{'numeric'},{'vector','real','finite','>=',0,'<=',1});
    assert(o.maxIterations>=1 && o.maxIterations==floor(o.maxIterations), ...
        'maxIterations must be a positive integer.');
end

function checkpoints=iteration_checkpoints(opts)
% Include zero and the final iteration so each run has endpoint states.
    checkpoints=unique([0 round(opts.checkpointFractions(:)'*opts.maxIterations) ...
        opts.maxIterations]);
end
