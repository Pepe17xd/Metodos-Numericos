function resultado = solve_thermoplate_cg(A, b, x0, opciones)
%SOLVE_THERMOPLATE_CG Compatibilidad; use metodo_gradiente_conjugado.
    resultado = metodo_gradiente_conjugado(A, b, x0, opciones);
end
