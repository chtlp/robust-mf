function [gU, gV] = ALM_grad_UV(M, W, U, V, T, s1, s2, a1, a2, mu)
    UV = U*V';
    gU = W .* (a1 - a2) * V + mu * W .* (s1 - s2 + 2*UV - 2*M) * V;
    gV = W' .* (a1 - a2)' * U + mu * W' .* (s1 - s2 + 2*UV - 2*M)' * U;
end