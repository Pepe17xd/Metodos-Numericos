%% Thermoplate - escenario base y solucionadores
% Modelo estacionario 2-D de conduction con dos fuentes localizadas.
% Supuestos ilustrativos de simulation, no mediciones de una laptop.
clear; clc;
addpath(fileparts(mfilename('fullpath')));
cfg = load_thermoplate_config();
model = build_thermoplate_model(cfg);
initial = cfg.initialTemperature*ones(model.nUnknowns,1);
opts = struct('tolerance',cfg.tolerance,'maxIterations',cfg.maxIterations);

%% Referencia directa dispersa
tic; x_direct = model.A\model.b; direct_time_s = toc;
direct_metrics = thermoplate_metrics(model,x_direct)

%% Iteration Gauss-Seidel y Gradiente Conjugado explicit
gs = solve_thermoplate_gs(model.A,model.b,initial,opts);
cg = solve_thermoplate_cg(model.A,model.b,initial,opts);
solver_summary = table(string({gs.solver;cg.solver}),[gs.iterations;cg.iterations], ...
    [gs.relativeResidual;cg.relativeResidual],[gs.time_s;cg.time_s], ...
    [string(gs.status);string(cg.status)], ...
    'VariableNames',{'solver','iter','rel_residual','time_s','status'})

%% Comprobaciones del caso base
assert(abs(model.sourcePower.cpu-cfg.cpu.power)<1e-10);
assert(abs(model.sourcePower.gpu-cfg.gpu.power)<1e-10);
assert(norm(gs.x-x_direct,inf)<2e-5 && norm(cg.x-x_direct,inf)<1e-5);
temperature_field = thermoplate_field(model,x_direct);
