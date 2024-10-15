% PRE: La cadena de Markov está implementada como un grafo de la siguiente forma:
% 'entrantes' es un cell array donde cada fila 'entrantes{i}' contiene una lista con los índices de las páginas que enlazan a la página i.
% 'salientes' es un array donde cada elemento 'salientes(j)' es un entero que indica la cantidad de enlaces salientes de la página j.
% 'alpha' es un número con 0 < alpha <= 1 que indica la probabilidad de que un usuario siga un enlace en la página actual. DEFAULT = 0.85
% 'max_iter' es un entero positivo que indica la cantidad máxima de iteraciones que se realizarán. DEFAULT = 100
% 'tol' es un número positivo que indica la tolerancia para el criterio de parada. DEFAULT = 1e-7

% POS: x_GS es el vector PageRank obtenido con el método de Gauss-Seidel y errors_GS es un vector que contiene el error en cada iteración.
function [x_GS, errors_GS] = gauss_seidel(entrantes, salientes, alpha, max_iter, tol)
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
    x_GS = ones(n, 1) / n;
    err_GS = tol + 1;
    iter_GS = 0;
    errors_GS = [];

    v = zeros(n, 1) + 1 / n;

    % Construir la matriz de transición

    alphaP = alpha * construirP(entrantes, salientes);;

    % A = [I - alphaP] => D = I
    D = sparse(eye(n));
    F = triu(alphaP, 1);
    E = tril(alphaP, -1);

    % Hallamos el vector PageRank con el método de Gauss-Seidel.
    while (err_GS > tol && iter_GS < max_iter)
        xant = x_GS;
        x_GS = (D - E) \ (F * x_GS + v);

        err_GS = norm(x_GS - xant);
        errors_GS = [errors_GS; err_GS];
        iter_GS += 1;
    end

    x_GS = x_GS / norm(x_GS, 1);
    toc
    err_GS
    iter_GS
endfunction
