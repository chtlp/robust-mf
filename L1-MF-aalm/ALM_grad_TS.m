function [gT, gs1, gs2] = ALM_grad_TS(M, W, U, V, T, s1, s2, a1, a2, mu)
    UV = U*V';
    L1 = s1 + UV - M - T;
    L2 = s2 + M - UV - T;

    gT = W .* (1 - a1 - a2 - mu * L1 - mu * L2);
    gs1 = W .* (a1 + mu * L1);
    gs2 = W .* (a2 + mu * L2);

    eps = 1e-6;
    gs1 = gs1 .* (~(s1 <= eps & gs1 > 0));
    gs2 = gs2 .* (~(s2 <= eps & gs2 > 0));
end