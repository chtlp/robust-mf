function L = lagrangian(M, W, U, V, S, Y, mu)

D = M - U*V' - S;
L = abs(S) + Y .* D + mu/2 * D.^2;
L = sum(sum(W .* L));

end