% update the primal variables
function [U1, V1, T, s1, s2] = ALM_update_primal(M, W, U, V, T, s1, s2, a1, a2, mu)

[m, n] = size(M);
k = size(U, 2);

UV = U*V';

[gT, gs1, gs2] = ALM_grad_TS(M, W, U, V, T, s1, s2, a1, a2, ...
                             mu);
g = max(max(abs(gT) + abs(gs1) + abs(gs2)));
L = lagrangian(M, W, U, V, T, s1, s2, a1, a2, mu);
fprintf('max(|gT| + |gs1| + |gs2|) = %.3f, L = %.3f\n', g, L);

for iter = 1:1
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

[gradU, gradV] = ALM_grad_UV(M, W, U, V, T, s1, s2, a1, a2, mu);
g = max(max(abs(gradU))) + max(max(abs(gradV)));
L = lagrangian(M, W, U, V, T, s1, s2, a1, a2, mu);

fprintf('max(|gU| + |gV|) = %.3f, L = %.3f\n', g, L);

for iter = 1:1
    % m x k
    Y = W .* (a1-a2+mu*(s1-s2-2*M)) * V / (-2*mu);
    % U1 = zeros(m, k);
    % for i = 1:m
    %     y = Y(i,:);
    %     w = W(i,:);
    %     R = V(w,:);
    %     A = R' * R;
    %     U1(i,:) = y / A;
    % end
    U1 = Y / (V' * V);

    % R = Y - W .* (U1 * V') * V;
    % assert(norm(vec(R), Inf) < 1e-2);

    % [gradU, gradV] = ALM_grad_UV(M, W, U1, V, T1, s1, s2, a1, a2, mu);
    % disp(gradU);
    % assert(norm(vec(gradU), Inf) < 1e-2);


    Y = W' .* (a1-a2+mu*(s1-s2-2*M))' * U1 / (-2*mu);
    % V1 = zeros(n, k);
    % for i = 1:n
    %     y = Y(i,:);
    %     w = W(:,i);
    %     R = U1(w,:);
    %     A = R' * R;
    %     V1(i,:) = y / A;
    % end
    V1 = Y / (U1' * U1);

    % [gradU, gradV] = ALM_grad_UV(M, W, U1, V1, T1, s1, s2, a1, a2, mu);
    % assert(norm(vec(gradV), Inf) < 1e-2);
    U = U1; V = V1;

    [gradU, gradV] = ALM_grad_UV(M, W, U, V, T, s1, s2, a1, a2, mu);
    g = max(max(abs(gradU))) + max(max(abs(gradV)));
    L = lagrangian(M, W, U, V, T, s1, s2, a1, a2, mu);

    fprintf('max(|gU| + |gV|) = %.3f, L = %.3f\n', g, L);
    if g < 1e-5
        break
    end

end

end
