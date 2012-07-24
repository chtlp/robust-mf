
U = rand(m, k);
V = rand(n, k);

a1 = zeros(m, n);
a2 = zeros(m, n);

T = abs(M - U*V');
s1 = T - (U*V' - M);
s2 = T - (M - U*V');

mu = 10;
eps_uv = 1e-4;

for iter = 1:2
    fprintf('iter = %d\n', iter);

    for jter = 1:20
        for kter = 1:4
            [V1, T1, s1, s2] = AALM_update_VST(M, W, U, V, T, s1, ...
                                               s2, a1, a2, mu);
            V = V1; T = T1;
            [gT, gs1, gs2] = ALM_grad_TS(M, W, U, V, T, s1, s2, a1, a2, ...
                                         mu);
            g = max(max(abs(gT) + abs(gs1) + abs(gs2)));
            fprintf('kiter = %d, max(|gT| + |gs1| + |gs2|) = %.3f\n', kter, g);
            if g < 1e-5
                break
            end
        end
        
        UV = U*V';
        a1 = a1 + mu * (s1 + UV - M - T);
        a2 = a2 + mu * (s2 + M - UV - T);        
        
        l1 = s1 + U*V' - M - T;
        l2 = s2 + M - U*V' - T;
        vio = sum(sum(abs(l1) + abs(l2)));
        fprintf('jter = %d, violation = %.3f\n', jter, vio);        
    end

    U = AALM_update_U(M, W, U, V, T, s1, s2, a1, a2, mu);
    
    fprintf('obj = %.3f\n', sum(sum(W .* abs(M - U*V'))));
    fprintf('\n');
end

