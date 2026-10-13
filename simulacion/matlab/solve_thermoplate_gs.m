function resultado = solve_thermoplate_gs(A, b, x0, opciones)
%SOLVE_THERMOPLATE_GS Compatibilidad; use metodo_gauss_seidel.
    resultado = metodo_gauss_seidel(A, b, x0, opciones);
end
