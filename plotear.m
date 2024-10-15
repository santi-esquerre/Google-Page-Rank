## -*- texinfo -*-
## @node plotear
## @deftypefn {Function} {@var{x} =} plotear ()
##
## @end deftypefn
function x = plotear ()

    clear all
    cd 'datos adjuntos'
    load -binary entrantes_full.dat
    load -binary salientes_full.dat
    cd ..

    % Calcular errores para cada método
    [~, errors_J] = jacobi(entrantes, salientes, 0.85, 150, 1e-16);
    [~, errors_GS] = gauss_seidel(entrantes, salientes, 0.85, 150, 1e-16);
    [~, errors_SOR] = SOR(entrantes, salientes, 0.85, 150, 1e-16);
    [~, errors_pot] = metodopotenciav2(entrantes, salientes, 0.85, 150, 1e-16);

    % Graficar los errores en función de las iteraciones
    figure;
    semilogy(1:length(errors_J), errors_J, '-', 'DisplayName', 'Jacobi');
    hold on;
    semilogy(1:length(errors_GS), errors_GS, '-', 'DisplayName', 'Gauss-Seidel');
    semilogy(1:length(errors_SOR), errors_SOR, '-', 'DisplayName', 'SOR');
    semilogy(1:length(errors_pot), errors_pot, '-', 'DisplayName', 'Potencia');
    xlabel('Iteraciones');
    ylabel('Error (norma)');
    title('Convergencia de los métodos iterativos');
    legend show;
    grid on;
endfunction
