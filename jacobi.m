% PRE: La cadena de Markov está implementada como un grafo de la siguiente forma:
% 'entrantes' es un cell array donde cada fila 'entrantes{i}' contiene una lista con los índices de las páginas que enlazan a la página i.
% 'salientes' es un array donde cada elemento 'salientes(j)' es un entero que indica la cantidad de enlaces salientes de la página j.
% 'alpha' es un número con 0 < alpha <= 1 que indica la probabilidad de que un usuario siga un enlace en la página actual. DEFAULT = 0.85
% 'max_iter' es un entero positivo que indica la cantidad máxima de iteraciones que se realizarán. DEFAULT = 100
% 'tol' es un número positivo que indica la tolerancia para el criterio de parada. DEFAULT = 1e-7

% POS: x_J es el vector PageRank obtenido con el método de las Jacobi y errors_J es un vector que contiene el error en cada iteración.
function [x_J, errors_J] = jacobi (entrantes, salientes, alpha, max_iter, tol)
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
    x_J = ones(n, 1) / n;
    err_J = tol + 1;
    iter_J = 0;
    errors_J = [];

    v = zeros(n, 1) + 1 / n;

    % Construir la matriz de transición
    alphaP = alpha * construirP(entrantes, salientes);

    % Aprovechamos que P tiene ceros en su diagonal, y por lo tanto D = inv(D) = I y ademas Q_J = inv(D)*(E+F) = I*alphaP = alphaP
    Q_J = alphaP;

    % r_J = inv(D)*v = I*v = v
    r_J = v;

    % Hallamos el vector PageRank con el método de Jacobi.
    while (err_J > tol && iter_J < max_iter)
        xant = x_J;
        x_J = Q_J * x_J + r_J;

        err_J = norm(x_J - xant);
        errors_J = [errors_J; err_J];
        iter_J += 1;
    end

    x_J = x_J / norm(x_J, 1);
    toc
    err_J
    iter_J
endfunction
