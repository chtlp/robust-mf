% update U
function [U1] = AALM_update_U(M, W, U, V, T, s1, s2, a1, a2, mu)
[m, k] = size(U);
    
Y = W .* (a1-a2+mu*(s1-s2-2*M)) * V / (-2*mu);
U1 = zeros(m, k);
for i = 1:m
    y = Y(i,:);
    w = W(i,:);
    A = V(w,:)' * V(w,:);
    U1(i,:) = y / A;
end

end
