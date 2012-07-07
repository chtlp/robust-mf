function U1 = update_U(M, W, U, V)

[m, n] = size(M);
k = size(U, 2);

U1 = zeros(m, k);

for i = 1:m
    w = W(i,:);
    y = M(i,w);
    A = V(w,:)';
    U1(i,:) = y / A;
end

end