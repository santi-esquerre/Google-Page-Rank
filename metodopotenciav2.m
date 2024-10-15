% PRE: La cadena de Markov está implementada como un grafo de la siguiente forma:
% 'entrantes' es un cell array donde cada fila 'entrantes{i}' contiene una lista con los índices de las páginas que enlazan a la página i.
% 'salientes' es un array donde cada elemento 'salientes(j)' es un entero que indica la cantidad de enlaces salientes de la página j.
% 'alpha' es un número con 0 < alpha <= 1 que indica la probabilidad de que un usuario siga un enlace en la página actual. DEFAULT = 0.85
% 'max_iter' es un entero positivo que indica la cantidad máxima de iteraciones que se realizarán. DEFAULT = 100
% 'tol' es un número positivo que indica la tolerancia para el criterio de parada. DEFAULT = 1e-7

% POS: x_POT es el vector PageRank obtenido con el método de las potencias y errors_POT es un vector que contiene el error en cada iteración.
function [x_POT, errors_POT] = metodopotenciav2 (entrantes, salientes, alpha, max_iter, tol)

    tic; % Controlar tiempo de ejecución

    % Valores por defecto
    if (nargin == 2)
        alpha = 0.85;
        max_iter = 100;
        tol = 1e-7;
    elseif (nargin == 3)
        max_iter = 100;
        tol = 1e-7;
    elseif (nargin == 4)
        tol = 1e-7;
    end

    n = length(salientes); % Número de páginas

    % Inicializamos variables para la iteración
    x_POT = ones(n, 1) / n;
    err_POT = tol + 1;
    iter_POT = 0;
    errors_POT = [];

    d = find(salientes == 0);
    vd = zeros(n, 1);
    vd(d) = 1 / n;

    % Construir matriz de transición
    P = construirP(entrantes, salientes);

    % Hallamos el vector PageRank con el método de las potencias.
    while (err_POT > tol && iter_POT < max_iter)
        xv = x_POT;

        % Cálculo vectorizado de Px
        Px = P * x_POT;

        % Cálculo vectorizado de vdx
        vdx = sum(vd .* x_POT);

        % Unimos los cálculos realizados
        x_POT = (alpha * (Px + vdx)) + ((1 - alpha) / n);

        err_POT = norm(x_POT - xv);
        errors_POT = [errors_POT; err_POT];
        ++iter_POT;
    end

    toc
    err_POT
    iter_POT
endfunction
