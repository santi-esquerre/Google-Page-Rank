% Cargamos datos
clear all
cd 'datos adjuntos'
load -binary entrantes_v1.dat
load -binary salientes_v1.dat
cd ..
entrantes
salientes

% Matriz de transicion original solo usando enlaces
n = size(salientes, 2); % Numero de paginas
P = zeros(n);

for i = 1:n
    enlaces = entrantes{i};
    P(i, enlaces) = 1 ./ salientes(enlaces);
endfor

P
% Agregamos links a paginas sin salida
d = find (salientes == 0);
P1 = P;
P1(:, d) = 1 / n

% Matriz de transicion final con factor de amortiguamiento
alpha = 0.85;
P2 = alpha * P1 + (1 - alpha) / n

% Hallamos vector PageRank con el metodo de las potencias
x = ones(n, 1) / n;
max_iter = 100;
tol = 1e-10;
error = tol + 1;

while error > tol
    xv = x;
    x = P2 * x;
    error = norm(x - xv);
endwhile

x
error
