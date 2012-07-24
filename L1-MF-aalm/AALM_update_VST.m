% update the primal variables
function [V1, T, s1, s2] = AALM_update_VST(M, W, U, V, T, s1, s2, a1, a2, mu)

[m, n] = size(M);
k = size(U, 2);

UV = U*V';

for iter = 1:2
    s1 = max(0, - (a1 + mu * (UV - M - T)) / mu);
    s2 = max(0, - (a2 + mu * (M - UV - T)) / mu);
    T = (1 - a1 - a2 - mu*s1 - mu*s2) / (-2*mu);

    [gT, gs1, gs2] = ALM_grad_TS(M, W, U, V, T, s1, s2, a1, a2, ...
                                 mu);
    g = max(max(abs(gT) + abs(gs1) + abs(gs2)));
    L = lagrangian(M, W, U, V, T, s1, s2, a1, a2, mu);
    fprintf('max(|gT| + |gs1| + |gs2|) = %.3f, L = %.3f\n', g, L);
    if g < 1e-5
        break
    end
end

U1 = U;

Y = W' .* (a1-a2+mu*(s1-s2-2*M))' * U1 / (-2*mu);
V1 = zeros(n, k);
for i = 1:n
    y = Y(i,:);
    w = W(:,i);
    R = U1(w,:);
    A = R' * R;
    V1(i,:) = y / A;
end

end
