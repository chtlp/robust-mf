
U = rand(m, k) * 0.1;
V = rand(n, k) * 0.1;

mu = 0.1;
for iter = 1:10
    [U1, mu] = gradL2_update_U(M, W, U, V, mu);
    [V1, mu] = gradL2_update_U(M', W', V, U1, mu);

    U = U1; V = V1;
    E = W .* (U*V' - M);
    err = sum(sum(E.^2));

    fprintf('iter = %d, obj = %.3f, mu = %.3e\n', iter, err, mu);
    if mu < 1e-7
        break
    else
        mu = mu * 1.05;
    end
end