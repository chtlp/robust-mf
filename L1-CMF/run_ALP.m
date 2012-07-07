opt = struct();

cvx_precision low

errs = norm(vec(E .* W), 1);
fprintf('model error %.3f\n', errs);

U = rand(m, k);
V = rand(n, k);
for iter=1:10
    err = norm(vec((U*V' - M) .* W), 1);
    fprintf('iter %d, training error %.3f\n', iter-1, err);
    U1 = comp_U_from_V(U, V, M, W, opt);
    fprintf('|U1-U| = %.3f\n', norm(vec(U1-U), 1));
    U = U1;
    
    V1 = comp_U_from_V(V, U, M', W', opt); 
    fprintf('|V1-V| = %.3f\n', norm(vec(V1-V), 1));
    V = V1;
end
