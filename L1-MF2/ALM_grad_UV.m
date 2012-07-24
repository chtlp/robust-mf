function [gU, gV] = ALM_grad_UV(M, W, U, V, S, Y, mu)
    UV = U*V';
    gU = - W .* Y * V - mu * W .* (M - UV - S) * V;
    gV = - W' .* Y' * U - mu * W' .* (M - UV - S)' * U;
end