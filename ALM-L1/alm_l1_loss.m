function L = alm_l1_loss(d, a1, a2, mu)

cvx_begin quiet

variable s1; 
variable s2; 
variable t;
minimize(t + a1*(d+s1-t)+ mu/2*(d+s1-t)^2 + a2*(-d+s2-t) + ...
         mu/2*(-d+s2-t)^2);
subject to 
s1>=0; s2>=0;

cvx_end

L = cvx_optval;
end