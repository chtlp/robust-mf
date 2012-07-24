function U1 = comp_U_from_V(U, V, M, W, opt)

[m, k] = size(U);
[n, k] = size(V);
assert( all([m,n] == size(M) &  [m, n] == size(W)) );

U1 = zeros(m, k);

nline=0;
for i = 1:m
    w = W(i,:);
    A = V(w,:);
    b = M(i,w)';
    x = comp_x_from_Ab(A, b);
    U1(i,:) = x;

    msg = sprintf('i = %d', i);
    fprintf(repmat('\b',1,nline));
    fprintf(msg);
    nline=numel(msg);
end
fprintf('\n');

end

function x = comp_x_from_Ab(A, b)

% disp(size(A));
% disp(size(b));

[m, n] = size(A);
assert( all([m, 1] == size(b)) );

cvx_begin quiet
    variable x(n);
    minimize( norm(A*x-b, 1) );
cvx_end

end