%NGDP_steady_state
%Insert used in NGDP_2024_Baseline.mod and its variants

coef_a = (betta*chi)^(1/(1-eps))*(1+n)^( eps/(1-eps) )*phi*(chi-1)^2/ybar^2;
coef_b = chi*(1 + (betta*chi)^(1/(1-eps))*(1+n)^( eps/(1-eps) )*(1 + 2*phi*(chi-1)*gbar/(chi*ybar^2)) );
coef_c = alfa*ybar + (betta*chi)^(1/(1-eps))*(1+n)^( eps/(1-eps) )*( (1 + phi*gbar/ybar^2 )*gbar - (1-alfa)*ybar);

p = [coef_a coef_b coef_c];
rts = roots(p);
b_root = max(rts);