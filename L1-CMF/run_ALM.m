
U = rand(m, k);
V = rand(n, k);

a1 = zeros(m, n);
a2 = zeros(m, n);

T = abs(M - U*V');
s1 = T - (U*V' - M);
s2 = T - (M - U*V');

mu = 10;
eps_uv = 1e-4;

for iter = 1:20
    fprintf('iter = %d\n', iter);


    for jter = 1:4
        [U1, V1, T1, s1, s2] = ALM_update_primal(M, W, U, V, T, s1, ...
                                                 s2, a1, a2, mu);
        U = U1; V = V1; T = T1;
        [gT, gs1, gs2] = ALM_grad_TS(M, W, U, V, T, s1, s2, a1, a2, ...
                                     mu);
        g = max(max(abs(gT) + abs(gs1) + abs(gs2)));
        fprintf('jiter = %d, max(|gT| + |gs1| + |gs2|) = %.3f\n', jter, g);
        if g < 1e-5
            break
        end
    end

    l1 = s1 + U*V' - M - T;
    l2 = s2 + M - U*V' - T;
    vio = sum(sum(abs(l1) + abs(l2)));
    fprintf('violation = %.3f\n', vio);
    fprintf('obj = %.3f\n', sum(sum(W .* abs(M - U*V'))));
    fprintf('\n');

    UV = U*V';
    a1 = a1 + mu * (s1 + UV - M - T);
    a2 = a2 + mu * (s2 + M - UV - T);

    % if iter == 10
    %     mu = mu * 5;
    % end
    
end

