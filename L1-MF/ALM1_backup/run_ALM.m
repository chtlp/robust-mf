
U = rand(m, k);
V = rand(n, k);

a1 = zeros(m, n);
a2 = zeros(m, n);

T = abs(M - U*V');
% s1 = max(0, - (a1 + mu * (U*V' - M - T)) / mu);
% s2 = max(0, - (a2 + mu * (M - U*V' - T)) / mu);
% T = (1 - a1 - a2 - mu*s1 - mu*s2) / (-2*mu);



mu = 500;
eps_uv = 1e-4;

for iter = 1:10
    fprintf('iter = %d\n', iter);

    % s1 = max(0, - (a1 + mu * (U*V' - M - T)) / mu);
    % s2 = max(0, - (a2 + mu * (M - U*V' - T)) / mu);
    % T = (1 - a1 - a2 - mu*s1 - mu*s2) / (-2*mu);

    for jter = 1:20
        [U1, V1, T1, s1, s2] = ALM_update_primal(M, W, U, V, T, s1, ...
                                                 s2, a1, a2, mu);
        [gradU, gradV] = ALM_grad_UV(M, W, U1, V1, T1, s1, s2, a1, a2, mu);
        U = U1; V = V1; T = T1;

        l1 = s1 + U*V' - M - T;
        l2 = s2 + M - U*V' - T;
        L = W .* (T + a1 .* l1 + a2 .* l2 + mu/2 * (l1 .* 2) + mu/2 * (l2 .* 2));
        ngU = norm(gradU, 'fro');
        ngV = norm(gradV, 'fro');
        fprintf('\tL = %.3f, fro(gU) = %.3f, fro(gV) = %.3f\n', ...
                sum(sum(L)), ngU, ngV);

        max_gu = max(max(abs(gradU)));
        max_gv = max(max(abs(gradV)));


        if max_gu < eps_uv && max_gv < eps_uv
            break
        end 
    end

    vio = sum(sum(abs(l1) + abs(l2)));
    fprintf('violation = %.3f\n', vio);
    fprintf('obj = %.3f\n', sum(sum(W .* abs(M - U*V'))));
    fprintf('\n');

    UV = U*V';
    a1 = a1 + mu * (s1 + UV - M - T);
    a2 = a2 + mu * (s2 + M - UV - T);

end

