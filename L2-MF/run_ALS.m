
U = rand(m, k);
V = rand(n, k);

for iter = 1:10
    U1 = update_U(M, W, U, V);
    V1 = update_U(M', W', V, U1);

    U = U1; V = V1;
    E = W .* (U*V' - M);
    err = sum(sum(E.^2));

    gU = E * V;
    gV = E' * U;
    fprintf('iter = %d, error = %.3f\n', iter, err);
    fprintf('\tgU = %.3f gV = %.3f\n', norm(gU, 'fro'), norm(gV, 'fro'));
end