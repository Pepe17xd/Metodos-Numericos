%% Thermoplate - execution reproducible de experimentos
% Escribe tablas CSV, metadatos JSON y estructura MAT en simulacion/resultados.
% Iteration runtime is computational cost, not physical heating time.
clear; clc;
addpath(fileparts(mfilename('fullpath')));
report = ejecutar_experimentos();
disp(report);
