m = 100;
n = 90;
k = 6;

% U-star, V-star, M-star
Us = rand(m, k);
Vs = rand(n, k);
Ms = Us*Vs';

% error matrix
E = std(vec(Ms)) * rand(m, n) .* (rand(m, n) < 0.1);

% 0.5 of the data available
W = rand(m, n) < 0.5;

M = Ms + E;

model_error = sum(sum(abs(E.*W)));


