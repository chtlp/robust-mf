% update the primal variables
function [U1, V1, S1] = ALM_update_primal(M, W, U, V, S, Y, mu, max_iter)

[m, n] = size(M);
k = size(U, 2);

S1 = M - U*V' + Y/mu;
S1 = sign(S1) .* max(0, abs(S1) - 1/mu);

for iter = 1:max_iter
    % m x k
    Z = M - S1 + Y ./ mu;
    U1 = zeros(m, k);
    for i = 1:m
        w = W(i,:);
        z = Z(i,w);
        A = V(w,:)';
        U1(i,:) = (z*A') / (A*A');
    end

    V1 = zeros(n, k);
    for i = 1:n
        w = W(:,i);
        z = Z(w,i)';
        A = U1(w,:)';
        V1(i,:) = (z*A') / (A*A');
    end

    U = U1; V = V1;
    [gradU, gradV] = ALM_grad_UV(M, W, U, V, S, Y, mu);
    g = max(max(abs(gradU))) + max(max(abs(gradV)));
    L = lagrangian(M, W, U, V, S, Y, mu);
    fprintf('max(|gU|) + max(|gV|) = %.3f, L = %.3f\n', g, L);
end

end
