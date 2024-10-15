% PRE: La cadena de Markov está implementada como un grafo de la siguiente forma:
% 'entrantes' es un cell array donde cada fila 'entrantes{i}' contiene una lista con los índices de las páginas que enlazan a la página i.
% 'salientes' es un array donde cada elemento 'salientes(j)' es un entero que indica la cantidad de enlaces salientes de la página j.
% POS: P es una matriz de transición de la cadena de Markov dada.
function P = construirP (entrantes, salientes)

    n = length(salientes); % Número de páginas

    % Calcular el número total de elementos no nulos
    nnz = sum(cellfun(@length, entrantes));

    % Preasignar vectores
    I = zeros(nnz, 1); % Indica las filas
    J = zeros(nnz, 1); % Indica las columnas
    V = zeros(nnz, 1); % Indica los valores (V(i) = P(I(i), J(i)))

    % Generar índices y valores
    idx = 1;

    for i = 1:n
        len = length(entrantes{i});

        if (len > 0)
            I(idx:idx + len - 1) = i;
            J(idx:idx + len - 1) = entrantes{i};
            V(idx:idx + len - 1) = 1 ./ salientes(entrantes{i});
            idx = idx + len;
        end

    end

    % Construir la matriz dispersa
    P = sparse(I, J, V, n, n);
endfunction
