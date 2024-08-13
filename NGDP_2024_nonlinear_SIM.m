%NGDP_2024_nonlinear_SIM
%Insert for the file NGDP_2024_nonlinear_LOOP.m
%Nonlinear simulation of the simple model (Algorithm in Supp Appendix).
%Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). Any errors are my own.

pi = NaN(T_sim,1); tau = pi; c1 = pi; c2 = pi; Index = pi; Max_Resid = pi;
R = pi; EV = pi; EV_reverse = pi; Utility = pi; ftau = pi; Resid_check = pi;
pi_prime = NaN(n_states^2,1); cprime = pi_prime; c_power = pi_prime; sdf_adj = pi_prime; 

for t=1:T_sim

    if t==1

        R_lag = R_init;
        y_lag = y_init;

    else

        R_lag = R(t-1);
        y_lag = y(t-1);

    end

    pi(t) = pistar*exp(e_vec(t))/(y(t)/y_lag)^(1-dummy_IT);
    tau(t) = ( gbar + (R_lag/((1+n)*pi(t)) - 1)*bstar ) / y(t);
    ftau(t) = phi*tau(t)^2;
    c1(t) = (1-alfa)*y(t) - tau(t)*y(t) - phi*tau(t)^2 - bstar;
    c2(t) = alfa*(1+n)*y(t) + R_lag/pi(t)*bstar;

    pi_prime = pistar*exp(states(:,1))*y(t)^(1-dummy_IT).*y_prime.^(-(1-dummy_IT));
    Resid0 = NaN(N_guess0,1);  Resid = NaN(N_guess,1);

    for k0=1:N_guess0

        Dum = 0;
        R_guess = R_guess_stack0(k0);    
        cprime0 = alfa*(1+n)*y_prime + R_guess*bstar.*pi_prime.^(-1);
        c_power0 = cprime0.^(1-gama);

        E_c_power0 = prob*c_power0;
        E_reverse0 = E_c_power0^(1/(1-gama));

        sdf_adj0 = betta*c1(t)^(1-eps).*cprime0.^(eps-1)*( 1 / E_reverse0 )^(1-gama-eps).*cprime0.^(1-gama-eps).*pi_prime.^(-1);

        Resid0(k0) = abs(1 - R_guess*prob*sdf_adj0);

        if k0 > 1 && Resid0(k0) > Resid0(k0-1)
            Dum = 1;
            break 
        end  

    end

        [Resid_mini0,Index_min0] = min(Resid0);
        R_guess_stack = R_guess_stack0(Index_min0) + linspace(R_l,R_u,N_guess);

    for k=1:N_guess 

        R_guess = R_guess_stack(k);
        cprime = alfa*(1+n)*y_prime + R_guess*bstar.*pi_prime.^(-1);
        c_power = cprime.^(1-gama);

        E_c_power = prob*c_power;
        E_reverse = E_c_power^(1/(1-gama));

        sdf_adj = betta*c1(t)^(1-eps).*cprime.^(eps-1)*( 1 / E_reverse )^(1-gama-eps).*cprime.^(1-gama-eps).*pi_prime.^(-1);

        Resid(k) = abs(1 - R_guess*prob*sdf_adj);

        if k>1 && Resid(k) > Resid(k-1)
            break 
        end  

    end

    [Resid_mini,Index_min] = min(Resid);
    R(t) = R_guess_stack(Index_min);

    Max_Resid(t) = Resid_mini;

    Index(t) = Index_min;
        
    cprime = alfa*(1+n).*y_prime + R(t)*bstar.*pi_prime.^(-1);
    c_power = cprime.^(1-gama);
    E_c_power = prob*c_power;
    E_reverse = E_c_power^(1/(1-gama));

    SDF_check = betta*c1(t)^(1-eps).*cprime.^(eps-1)*( 1 / E_reverse )^(1-gama-eps).*cprime.^(1-gama-eps).*pi_prime.^(-1);
    Resid_check(t) = abs(1 - R(t)*prob*SDF_check);

    EV(t) = prob*c_power;

    Utility(t) =  1/(1-gama)*( c1(t)^eps + betta*( EV(t) )^(eps/(1-gama)) )^((1-gama)/eps);

end




