function report = experimentos(outputRoot, timingOptions)
%EXPERIMENTOS Ejecuta caso base, refinamiento, sensibilidad y estres.
%   Optional TIMINGOPTIONS fields are warmupRuns and measuredRuns. Each
%   measured run starts from the same initial vector; its solution is not
%   used to replace the canonical solution/history calculated below.
%   Ejemplo: addpath('simulacion/matlab'); experimentos
    if nargin<1 || isempty(outputRoot)
        outputRoot=fullfile(fileparts(mfilename('fullpath')),'..','resultados');
    end
    if nargin<2 || isempty(timingOptions), timingOptions=struct(); end
    timingOptions=validate_timing_options(timingOptions);
    cfg=planteamiento_modelo();
    baseDir=fullfile(outputRoot,'escenario_base'); convDir=fullfile(outputRoot,'convergencia');
    refDir=fullfile(outputRoot,'refinamiento');
    senDir=fullfile(outputRoot,'sensibilidad'); stressDir=fullfile(outputRoot,'estres_numerico');
    dirs={baseDir,convDir,refDir,senDir,stressDir}; for i=1:numel(dirs),if ~isfolder(dirs{i}),mkdir(dirs{i});end,end

    % Scenario 1: same A, b and constant initial guess for all solvers.
    nx=cfg.baselineGrid(1); ny=cfg.baselineGrid(2); model=discretizacion_matriz(cfg,nx,ny);
    x0=cfg.initialTemperature*ones(model.nUnknowns,1);
    directTimer=tic; xDirect=model.A\model.b; directTime=toc(directTimer);
    opts=struct('tolerance',cfg.tolerance,'maxIterations',cfg.maxIterations, ...
        'checkpointFractions',cfg.checkpointFractions);
    gs=metodo_gauss_seidel(model.A,model.b,x0,opts);
    cg=metodo_gradiente_conjugado(model.A,model.b,x0,opts);
    [baseTiming,baseTimingStats]=benchmark_solvers(model.A,model.b,x0,opts, ...
        nx,ny,"base",timingOptions);
    [directTiming,directTimingStats]=benchmark_direct(model.A,model.b,nx,ny, ...
        "base",timingOptions);
    baseTiming=[directTiming;baseTiming];
    writetable(baseTiming,fullfile(baseDir,'timing_repeats.csv'));
    mDirect=metricas_temperatura(model,xDirect);
    field=campo_temperatura(model,xDirect); writetable(field,fullfile(baseDir,'baseline_field.csv'));
    history=history_table({gs,cg}); writetable(history,fullfile(baseDir,'convergence_history.csv'));
    exportar_estados_representativos(model,{gs,cg}, ...
        fullfile(baseDir,'representative_states.csv'));
    mGs=metricas_temperatura(model,gs.x); mCg=metricas_temperatura(model,cg.x);
    summary=summary_table(model,{gs,cg},[mGs.Tmax_C mCg.Tmax_C],baseTimingStats);
    directRow=table("Sparse direct",cfg.tolerance,1,mDirect.relative_residual,NaN, ...
        directTimingStats.median_s,directTime,directTimingStats.median_s, ...
        directTimingStats.dispersion_s,directTimingStats.iqr_s, ...
        directTimingStats.measured_runs, ...
        nx,ny, ...
        model.nUnknowns,cfg.cpu.power,cfg.gpu.power,mDirect.Tmax_C,"reference", ...
        'VariableNames',summary.Properties.VariableNames);
    summary=[directRow;summary]; writetable(summary,fullfile(baseDir,'solver_summary.csv'));
    sourceInfo=source_manifest(cfg);
    metadata=struct('experiment_id','scenario_base','model','steady_2d_conduction_localized_sources', ...
        'geometry',struct('Lx_m',cfg.Lx,'Ly_m',cfg.Ly,'sources',sourceInfo), ...
        'thickness_m',cfg.thickness,'k_W_mK',cfg.k,'boundary_temperature_C',cfg.Tb, ...
        'cpu_power_W',cfg.cpu.power,'gpu_power_W',cfg.gpu.power, ...
        'cpu_q_W_m3',cfg.cpu.power/(cfg.cpu.width*cfg.cpu.height*cfg.thickness), ...
        'gpu_q_W_m3',cfg.gpu.power/(cfg.gpu.width*cfg.gpu.height*cfg.thickness), ...
        'materialized_source_power_W',model.sourcePower,'nx',nx,'ny',ny, ...
        'unknowns',model.nUnknowns,'dx_m',model.dx,'dy_m',model.dy, ...
        'tolerance',cfg.tolerance,'assumptions',{cfg.assumptions}, ...
        'initial_state_C',cfg.initialTemperature,'matrix_symmetry_relative_error',model.symmetryError, ...
        'timing_method',struct('warmup_runs',timingOptions.warmupRuns, ...
        'measured_runs',timingOptions.measuredRuns, ...
        'dispersion','sample standard deviation of solver-reported elapsed seconds'), ...
        'matrix_spd_evidence','Verified by symmetry and successful sparse Cholesky factorization; Dirichlet energy is positive for nonzero interior vectors.');
    escribir_json(fullfile(baseDir,'metadata.json'),metadata);
    save(fullfile(baseDir,'baseline.mat'),'model','xDirect','gs','cg','metadata','-v7.3');
    write_manifest(baseDir,metadata,struct('field','baseline_field.csv', ...
        'history','convergence_history.csv','states','representative_states.csv', ...
        'mesh','solver_summary.csv','summary','solver_summary.csv', ...
        'timing','timing_repeats.csv','metadata','metadata.json','mat','baseline.mat'));
    % Dedicated convergence folder mirrors the iterative datasets while the
    % field and summary remain shared references to the base experiment.
    copyfile(fullfile(baseDir,'convergence_history.csv'),fullfile(convDir,'convergence_history.csv'),'f');
    copyfile(fullfile(baseDir,'representative_states.csv'),fullfile(convDir,'representative_states.csv'),'f');
    convMeta=metadata; convMeta.experiment_id='iterative_convergence';
    convMeta.initial_state_C=cfg.initialTemperature;
    escribir_json(fullfile(convDir,'metadata.json'),convMeta);
    write_manifest(convDir,convMeta,struct('history','convergence_history.csv', ...
        'states','representative_states.csv','metadata','metadata.json'));

    % Scenario 3: nested grids. Direct sparse reference at each level separates
    % algebraic error (iterative versus direct) from spatial refinement.
    refinement=table(); refinementTiming=table();
    refSol=cell(size(cfg.grids,1),1); refModels=cell(size(cfg.grids,1),1);
    for g=1:size(cfg.grids,1)
        gx=cfg.grids(g,1); gy=cfg.grids(g,2); mm=discretizacion_matriz(cfg,gx,gy);
        xd=mm.A\mm.b; md=metricas_temperatura(mm,xd); refSol{g}=reshape_field(mm,xd); refModels{g}=mm;
        ini=cfg.initialTemperature*ones(mm.nUnknowns,1);
        gg=metodo_gauss_seidel(mm.A,mm.b,ini,opts);
        cc=metodo_gradiente_conjugado(mm.A,mm.b,ini,opts);
        [meshTiming,meshTimingStats]=benchmark_solvers(mm.A,mm.b,ini,opts, ...
            gx,gy,"refinement",timingOptions);
        refinementTiming=[refinementTiming;meshTiming]; %#ok<AGROW>
        mg=metricas_temperatura(mm,gg.x); mc=metricas_temperatura(mm,cc.x);
        rows=summary_table(mm,{gg,cc},[mg.Tmax_C mc.Tmax_C],meshTimingStats);
        rows.spatial_max_abs_difference_C=NaN(height(rows),1);
        rows.algebraic_max_abs_error_C=[max(abs(gg.x-xd));max(abs(cc.x-xd))];
        rows.discrete_source_power_W=repmat(mm.sourcePower.cpu+mm.sourcePower.gpu,height(rows),1);
        refinement=[refinement;rows]; %#ok<AGROW>
    end
    for g=1:size(cfg.grids,1)-1
        coarse=refSol{g}; fine=refSol{g+1}; cm=refModels{g}; fm=refModels{g+1};
        [Xc,Yc]=meshgrid(cm.x,cm.y); [Xf,Yf]=meshgrid(fm.x,fm.y);
        fineAtCoarse=interp2(Xf,Yf,fine,Xc,Yc,'linear');
        delta=max(abs(coarse-fineAtCoarse),[],'all');
        idx=(g-1)*2+(1:2); refinement.spatial_max_abs_difference_C(idx)=delta;
    end
    writetable(refinement,fullfile(refDir,'refinement.csv'));
    writetable(refinementTiming,fullfile(refDir,'timing_repeats.csv'));
    refMeta=struct('experiment_id','mesh_refinement','comparison','max nodal absolute difference: each direct coarse field versus bilinear interpolation of next finer direct field', ...
        'geometry',metadata.geometry,'tolerance',cfg.tolerance);
    escribir_json(fullfile(refDir,'metadata.json'),refMeta);
    refMeta.timing_method=metadata.timing_method;
    write_manifest(refDir,refMeta,struct('mesh','refinement.csv', ...
        'timing','timing_repeats.csv','metadata','metadata.json'));

    % Scenario 4: vary each source independently; direct solve isolates the
    % physical model response from iterative stopping error.
    scales=[0 0.5 1 1.5]; sensitivity=table();
    for s=1:numel(scales)
        for sourceIdx=1:2
            variant=cfg; sourceName={'cpu','gpu'}; sourceName=sourceName{sourceIdx};
            variant.(sourceName).power=cfg.(sourceName).power*scales(s);
            sm=discretizacion_matriz(variant,nx,ny); sx=sm.A\sm.b; mt=metricas_temperatura(sm,sx);
            sensitivity=[sensitivity;table(string(sourceName),scales(s),variant.cpu.power,variant.gpu.power, ...
                mt.Tmax_C,mt.hotspot_x_m,mt.hotspot_y_m,mt.relative_residual, ...
                'VariableNames',{'varied_source','scale','cpu_power_W','gpu_power_W','Tmax_C', ...
                'hotspot_x_m','hotspot_y_m','rel_residual'})]; %#ok<AGROW>
        end
    end
    writetable(sensitivity,fullfile(senDir,'sensitivity.csv'));
    senMeta=struct('experiment_id','source_sensitivity','geometry',metadata.geometry, ...
        'scales',[scales],'procedure','Vary CPU or GPU power independently; hold the other source and all other parameters fixed.');
    escribir_json(fullfile(senDir,'metadata.json'),senMeta);
    write_manifest(senDir,senMeta,struct('sensitivity','sensitivity.csv','metadata','metadata.json'));

    % Scenario 5: tolerance sweep plus explicitly labelled extreme and
    % nonphysical source cases. These are numerical stress cases, not laptop data.
    stress=table(); stressTiming=table(); tolList=[1e-6 1e-8 1e-10 1e-12];
    for tol=tolList
        so=opts; so.tolerance=tol; so.maxIterations=max(cfg.maxIterations,100000);
        [repeatRows,repeatStats]=benchmark_solvers(model.A,model.b,x0,so, ...
            nx,ny,"baseline_tolerance",timingOptions);
        stressTiming=[stressTiming;repeatRows]; %#ok<AGROW>
        for alg=1:2
            if alg==1, rr=metodo_gauss_seidel(model.A,model.b,x0,so); else,rr=metodo_gradiente_conjugado(model.A,model.b,x0,so);end
            stress=[stress;stress_row("baseline_tolerance",tol,rr,repeatStats(alg))]; %#ok<AGROW>
        end
    end
    cases={"extreme_positive_1000x",1000,1000;"nonphysical_negative_cpu",-1,1};
    for c=1:size(cases,1)
        variant=cfg; variant.cpu.power=cfg.cpu.power*cases{c,2}; variant.gpu.power=cfg.gpu.power*cases{c,3};
        sm=discretizacion_matriz(variant,nx,ny); xx0=variant.initialTemperature*ones(sm.nUnknowns,1);
        so=opts; so.tolerance=1e-8;
        [repeatRows,repeatStats]=benchmark_solvers(sm.A,sm.b,xx0,so, ...
            nx,ny,string(cases{c,1}),timingOptions);
        stressTiming=[stressTiming;repeatRows]; %#ok<AGROW>
        for alg=1:2
            if alg==1,rr=metodo_gauss_seidel(sm.A,sm.b,xx0,so);else,rr=metodo_gradiente_conjugado(sm.A,sm.b,xx0,so);end
            stress=[stress;stress_row(cases{c,1},so.tolerance,rr,repeatStats(alg))]; %#ok<AGROW>
        end
    end
    writetable(stress,fullfile(stressDir,'stress.csv'));
    writetable(stressTiming,fullfile(stressDir,'timing_repeats.csv'));
    stressMeta=struct('experiment_id','numerical_stress','geometry',metadata.geometry, ...
        'timing_method',metadata.timing_method, ...
        'warning','Extreme and negative source powers are deliberately nonphysical numerical probes, not laptop estimates.');
    escribir_json(fullfile(stressDir,'metadata.json'),stressMeta);
    write_manifest(stressDir,stressMeta,struct('stress','stress.csv', ...
        'timing','timing_repeats.csv','metadata','metadata.json'));
    report=struct('baseDirectory',baseDir,'refinement',refinement,'sensitivity',sensitivity,'stress',stress);
end

function H=history_table(results)
    H=table(); for a=1:numel(results)
        r=results{a}; H=[H;table((0:numel(r.history)-1)',repmat(string(r.solver),numel(r.history),1), ...
            r.history(:),'VariableNames',{'iter','solver','rel_residual'})]; %#ok<AGROW>
    end
end
function S=summary_table(m,results,tmax,timingStats)
    S=table(); for j=1:numel(results)
        r=results{j}; timing=timingStats(j);
        S=[S;table(string(r.solver),r.tolerance,r.iterations,r.relativeResidual, ...
            r.relativeStepError, ...
            timing.median_s,r.time_s,timing.median_s,timing.dispersion_s, ...
            timing.iqr_s,timing.measured_runs,m.nx,m.ny,m.nUnknowns,m.cfg.cpu.power, ...
            m.cfg.gpu.power,tmax(j),string(r.status), ...
            'VariableNames',{'solver','tolerance','iter','rel_residual','rel_step_error','time_s', ...
            'canonical_time_s','time_median_s','time_std_s','time_iqr_s','timing_repeats', ...
            'nx','ny','unknowns','cpu_power_W','gpu_power_W','Tmax_C','status'})]; %#ok<AGROW>
    end
end

function options=validate_timing_options(options)
% One warm-up and three measured runs are the default timing protocol.
    defaults=struct('warmupRuns',1,'measuredRuns',3);
    names=fieldnames(defaults);
    for k=1:numel(names)
        if ~isfield(options,names{k}), options.(names{k})=defaults.(names{k}); end
        validateattributes(options.(names{k}),{'numeric'}, ...
            {'scalar','real','finite','integer','>=',0});
    end
    assert(options.measuredRuns>=3, ...
        'Use at least three measured repetitions for timing summaries.');
end

function [rows,stats]=benchmark_solvers(A,b,x0,opts,nx,ny,experiment,options)
% Benchmark each explicit iterative algorithm without retaining its output.
    labels={"Gauss-Seidel","Conjugate Gradient"};
    rows=table(); stats=repmat(struct('median_s',NaN,'dispersion_s',NaN, ...
        'iqr_s',NaN, ...
        'measured_runs',options.measuredRuns),numel(labels),1);
    for alg=1:numel(labels)
        samples=zeros(options.measuredRuns,1);
        for run=1:options.warmupRuns
            execute_solver(alg,A,b,x0,opts); %#ok<NASGU>
        end
        for run=1:options.measuredRuns
            result=execute_solver(alg,A,b,x0,opts);
            samples(run)=result.time_s;
            rows=[rows;table(string(experiment),labels{alg},nx,ny,numel(b), ...
                "measured",run,result.time_s,result.iterations, ...
                result.relativeResidual,string(result.status), ...
                'VariableNames',{'experiment','solver','nx','ny','unknowns', ...
                'run_type','repetition','time_s','iterations','relative_residual','status'})]; %#ok<AGROW>
        end
        stats(alg).median_s=median(samples);
        stats(alg).dispersion_s=std(samples,0);
        stats(alg).iqr_s=iqr(samples);
    end
end

function [rows,stats]=benchmark_direct(A,b,nx,ny,experiment,options)
% Repeat the sparse backslash reference solve using the same timing protocol.
    samples=zeros(options.measuredRuns,1);
    rows=table();
    for run=1:options.warmupRuns
        warmupSolution=A\b; %#ok<NASGU>
    end
    for run=1:options.measuredRuns
        timer=tic; x=A\b; elapsed=toc(timer);
        samples(run)=elapsed;
        residual=norm(b-A*x)/max(norm(b),double(norm(b)==0));
        rows=[rows;table(string(experiment),"Sparse direct",nx,ny,numel(b), ...
            "measured",run,elapsed,1,residual,"reference", ...
            'VariableNames',{'experiment','solver','nx','ny','unknowns', ...
            'run_type','repetition','time_s','iterations','relative_residual','status'})]; %#ok<AGROW>
    end
    stats=struct('median_s',median(samples),'dispersion_s',std(samples,0), ...
        'iqr_s',iqr(samples), ...
        'measured_runs',options.measuredRuns);
end

function result=execute_solver(algorithm,A,b,x0,opts)
% Keep the benchmark dispatch explicit so students can map it to each method.
    if algorithm==1
        result=metodo_gauss_seidel(A,b,x0,opts);
    else
        result=metodo_gradiente_conjugado(A,b,x0,opts);
    end
end
function T=reshape_field(m,x)
    T=m.cfg.Tb*ones(m.ny,m.nx); T(2:end-1,2:end-1)=reshape(x,m.ny-2,m.nx-2);
end
function src=source_manifest(cfg)
    src=struct('name',{'CPU','GPU'},'x_m',{cfg.cpu.x,cfg.gpu.x}, ...
        'y_m',{cfg.cpu.y,cfg.gpu.y},'width_m',{cfg.cpu.width,cfg.gpu.width}, ...
        'height_m',{cfg.cpu.height,cfg.gpu.height});
end
function write_manifest(folder,metadata,datasets)
    names=fieldnames(datasets); entries=struct();
    for k=1:numel(names), entries.(names{k})=datasets.(names{k}); end
    manifest=struct('experiment_id',metadata.experiment_id,'metadata',metadata,'datasets',entries);
    escribir_json(fullfile(folder,'manifest.json'),manifest);
end
function r=stress_row(caseName,tol,s,timing)
    r=table(string(caseName),tol,string(s.solver),string(s.status),s.iterations, ...
        s.relativeResidual,timing.median_s,s.time_s,timing.median_s, ...
        timing.dispersion_s,timing.iqr_s,timing.measured_runs,'VariableNames', ...
        {'case','tol','solver','status','iterations','residual','time_s', ...
        'canonical_time_s','time_median_s','time_std_s','time_iqr_s','timing_repeats'});
end

