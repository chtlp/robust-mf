m = 500;
n = 450;
k = 4;

% U-star, V-star, M-star
Us = rand(m, k);
Vs = rand(n, k);
Ms = Us*Vs';

% error matrix
Es = std(vec(Ms)) * randn(m, n) .* (rand(m, n) < 0.1);

% 0.5 of the data available
W = rand(m, n) < 0.5;

M = Ms + Es;

E = Es .* W;
model_error = sum(sum(E.^2));


