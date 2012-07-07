function [U1, mu] = subgrad_update_U(M, W, U, V, mu)

E = W .* (U*V' - M);
err = sum(sum(abs(E)));
gradU = sign(E) * V;
g = norm(gradU, 'fro')^2;


while mu >= 1e-7
    E1 = W .* ((U - mu*gradU) * V' - M);
    err1 = sum(sum(abs(E1)));

    if err < err1
        mu = mu * 0.5;
    else
        break;
    end
end

U1 = U - mu * gradU;

end