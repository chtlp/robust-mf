function [U1, V1, mu] = gradL2_updateUV(M, W, U, V, mu)

E = W .* (U*V' - M);
err = sum(sum(E.^2));
gradU = 2 * E * V;
gradV = 2 * E' * U;


while mu >= 1e-7
    E1 = W .* ((U - mu*gradU) * (V - mu*gradV)' - M);
    err1 = sum(sum(E1.^2));

    if err < err1
        mu = mu * 0.5;
    else
        break;
    end
end

U1 = U - mu * gradU;
V1 = V - mu * gradV;

end