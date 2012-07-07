% fixing V, we only try to update U
U = rand(m, k);
V = Vs;

T = zeros(m, n);
alpha1 = zeros(m, n);
alpha2 = zeros(m, n);

% learning rate
eta = 0.1;

for iter=1:100
    UV = U*V';
    alpha1 = max(0, alpha1 + W .* (UV-M-T) * eta);
    alpha2 = max(0, alpha2 + W .* (M-UV-T) * eta);
    
    gV = W' .* (alpha1'- alpha2') * U;
    V = V - eta * gV;

    gU = W .* (alpha1-alpha2) * V;
    U = U - eta * gU;

    T = max(0, T - eta * W .* (1 - alpha1 - alpha2));
    
    uv = U*V';
    violation = norm(vec(uv-M-T), 1);
    t = norm(vec(T), 1);
    obj = norm(vec(uv-M), 1);
    fprintf('iter = %.3f, violation = %.3f, |T| = %.3f, obj = %.3f\n', ...
            iter, violation, t, obj);
end