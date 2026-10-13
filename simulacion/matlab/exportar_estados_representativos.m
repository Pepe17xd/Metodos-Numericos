function statesTable = exportar_estados_representativos(model, canonicalResults, outputPath)
%EXPORTAR_ESTADOS_REPRESENTATIVOS Exporta estados inicial, medio y final.
%   Repeats each explicit solver from the same uniform initial condition,
%   with its canonical iteration count as the limit. The midpoint checkpoint
%   is therefore relative to the actual run length, not the configured limit.
%   The replay must converge to the canonical solution before any field is
%   exported. OUTPUTPATH is optional; when supplied, the table is written CSV.
%
%   The iteration index is an algebraic solver count. It is not physical time.

    if nargin < 3
        outputPath = '';
    end
    validateattributes(canonicalResults, {'cell'}, {'nonempty'});
    assert(isfield(model, 'cfg') && isfield(model.cfg, 'initialTemperature'), ...
        'El modelo debe incluir la temperatura inicial en model.cfg.');

    statesTable = table();
    [X, Y] = meshgrid(model.x, model.y);
    numberOfNodes = numel(X);
    initialVector = model.cfg.initialTemperature * ones(model.nUnknowns, 1);

    for index = 1:numel(canonicalResults)
        canonical = canonicalResults{index};
        options = struct('tolerance', canonical.tolerance, ...
            'maxIterations', canonical.iterations, 'checkpointFractions', 0.5);

        switch canonical.solver
            case 'Gauss-Seidel'
                replay = metodo_gauss_seidel(model.A, model.b, initialVector, options);
            case 'Conjugate Gradient'
                replay = metodo_gradiente_conjugado(model.A, model.b, initialVector, options);
            otherwise
                error('THERMOPLATE:UnknownSolver', ...
                    'No existe una repetición de estados para el método "%s".', canonical.solver);
        end

        assert(strcmp(replay.status, canonical.status) && ...
            replay.iterations == canonical.iterations, ...
            'La repetición de %s no coincide con el estado canónico.', canonical.solver);
        solutionScale = max(1, norm(canonical.x, inf));
        assert(norm(replay.x - canonical.x, inf) <= 100 * eps(solutionScale), ...
            'La solución repetida de %s difiere de la solución canónica.', canonical.solver);
        assert(numel(replay.states) == 3 && replay.states(1).iter == 0 && ...
            replay.states(2).iter > 0 && replay.states(3).iter == canonical.iterations, ...
            'La repetición de %s no produjo los estados inicial, medio y final.', canonical.solver);

        for stateIndex = 1:numel(replay.states)
            solverState = replay.states(stateIndex);
            temperature = model.cfg.Tb * ones(model.ny, model.nx);
            temperature(2:end-1, 2:end-1) = ...
                reshape(solverState.x, model.ny-2, model.nx-2);
            rows = table(repmat(solverState.iter, numberOfNodes, 1), ...
                repmat(string(canonical.solver), numberOfNodes, 1), ...
                X(:), Y(:), temperature(:), ...
                'VariableNames', {'iter', 'solver', 'x_m', 'y_m', 'T_C'});
            statesTable = [statesTable; rows]; %#ok<AGROW>
        end
    end

    if ~isempty(outputPath)
        writetable(statesTable, outputPath);
    end
end
