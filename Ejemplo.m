% Definimos la matriz P
P = [0.8, 0.1, 0.2; 0.1, 0.7, 0.2; 0.1, 0.2, 0.6]
% V es una matriz cuyas columnas son vectores propios de P,
% D la matriz diagonal semejante a P,
% V(:, i) es el vector propio asociado al valor propio D(i,i)
[V, D] = eig(P)
% x es un vector propio de P asociado al valor propio 1,
% con componentes no negativas, y que cumple que norm(x, 1) = 1
x = abs(V(:, 1)) / norm(V(:, 1), 1)
