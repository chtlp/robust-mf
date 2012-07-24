
U = rand(m, k);
V = rand(n, k);
S = zeros(m, n);
Y = zeros(m, n);

% mu = 3;
eps_uv = 1e-4;
uv_iter=1;

vio_seq = [];
obj_seq = [];
for iter = 1:25
    fprintf('iter = %d\n', iter);

    [U, V, S] = ALM_update_primal(M, W, U, V, S, Y, mu, uv_iter);
    D = M - U*V' - S;
    Y = Y + mu * D;
 
    vio = sum(sum(abs(D)));
    obj = sum(sum(W .* abs(M - U*V')));
    fprintf('violation = %.3f\n', vio);
    fprintf('obj = %.3f\n', obj);
    fprintf('\n');

    vio_seq = [vio_seq vio];
    obj_seq = [obj_seq obj];
end

