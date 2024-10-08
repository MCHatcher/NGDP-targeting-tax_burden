%Steady_state_LS2_insert
%Computes approximation to lump-sum taxes steady state (currently unused)

n_guess = 9E6;
tol1 = 9e-6;
l_min = 0.1; l_max = 0.9;
tau_min = 0.01; tau_max = 0.4;
Loss1 = NaN(n_guess,1); Dum1 = 0;

for i = 1:n_guess

        lab_guess = l_min + rand*(l_max-l_min);
        tau_guess = tau_min + rand*(tau_max-tau_min);
        
        w1 = (1-alfa)*(alfa/chi1)^(alfa/(1-alfa));
        TR1 = tau_p*(1+n)*w1*lab_guess;

        s1 = ( ( betta*chi1*(1/(1-lab_guess))^( (1-thetta)*eps ) )^(1/(1-thetta*eps))*(w1*lab_guess - tau_guess) + psi*tau_guess - TR1) / ( chi1 + ( betta*chi1*(1/(1-lab_guess))^( (1-thetta)*eps ) )^(1/(1-thetta*eps)) );
        k_tild1 = (alfa/chi1)^(1/(1-alfa))*lab_guess;
        b1 = s1 - (1+n)*k_tild1; 

        c11 = w1*lab_guess - s1 - tau_guess;
        c21 = chi1*s1 - psi*tau_guess + TR1;

        Resid11 = (1-thetta)/thetta*c11 / w1 - (1-lab_guess);
        Resid21 = (1+psi/(1+n))*tau_guess - gbar - (chi1/(1+n) -1)*b1 - TR1/(1+n); 

        Loss1(i) = Resid11^2 + Resid21^2; 

        if Loss1(i) < tol1
            Dum1 = 1;
        end

        if Dum1==1
            break
        end

end

b1_ss = b1; 
k1_ss = (1+n)*k_tild1; 
s1_ss = s1; 
tau1_ss = tau_guess;
c11_ss = c11; c21_ss = c21;
l1_ss = lab_guess;
utility1_ss = (1/(1-gama))*(  (c11_ss^thetta*(1-l_ss)^(1-thetta))^eps + betta*(c21_ss^thetta)^eps )^((1-gama)/eps);

          