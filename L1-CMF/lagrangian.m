function L = lagrangian(M, W, U, V, T, s1, s2, a1, a2, mu)

l1 = s1 + U*V' - M - T;
l2 = s2 + M - U*V' - T;
L = W .* (T + a1 .* l1 + a2 .* l2 + mu/2 * (l1.^22) + mu/2 * (l2.^2));

L = sum(sum(L));

end