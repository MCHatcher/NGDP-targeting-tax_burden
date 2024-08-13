%NGDP_2024_nonlinear_sim_transition 
%Insert for the file NGDP_2024_nonlinear_transition.m and its variants
%Nonlinear simulation of transition to NIT in the simple model (Algorithm in Supp Appendix).
%Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). Any errors are my own.

pi = NaN(T_fin,1); tau = pi; c1 = pi; c2 = pi; Index = pi; Max_Resid = pi; Lambda_NIT = pi;
pi_prime = NaN(n_states^2,1); cprime = pi_prime; c_power = pi_prime; sdf_adj = pi_prime; 
R = NaN(T_fin,1); EV = R; EV_reverse = R; Utility = R; U_IT = Utility; U_sum = R; 

Dummy_init = 1;

for t=1:T_sim

    dummy_IT = 1;

    if t==1

        R_lag = R_init;
        y_lag = y_init;

    else

        R_lag = R(t-1);
        y_lag = y(t-1);

    end

    pi(t) = pistar*exp(e_vec(t))/(y(t)/y_lag)^(1-dummy_IT);
    tau(t) = ( gbar + (R_lag/((1+n)*pi(t)) - 1)*bstar ) / y(t);
    c1(t) = (1-alfa)*y(t) - tau(t)*y(t) - phi*tau(t)^2 - bstar;
    c2(t) = alfa*(1+n)*y(t) + R_lag/pi(t)*bstar;  

    pi_prime = pistar*exp(states(:,1))*y(t)^(1-dummy_IT).*y_prime.^(-(1-dummy_IT));
    Resid0 = NaN(N_guess0,1);  Resid = NaN(N_guess,1); Resid_init = NaN(N_guess_init,1); 

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

    EV(t) = prob*c_power;

    Utility(t) =  1/(1-gama)*( c1(t)^eps + betta*( EV(t) )^(eps/(1-gama)) )^((1-gama)/eps);

end

    if dummy_NIT==1
        dummy_IT = 0;
    end

    if Announced == 1

    pi_prime = pistar*exp(states(:,1))*y(T_sim)^(1-dummy_IT).*y_prime.^(-(1-dummy_IT));

    for k_init=1:N_guess_init

        R_guess = R_guess_stack_init(k_init);  
        cprime_init = alfa*(1+n)*y_prime + R_guess*bstar.*pi_prime.^(-1);
        c_power_init = cprime_init.^(1-gama);

        E_c_power_init = prob*c_power_init;
        E_reverse_init = E_c_power_init^(1/(1-gama));

        sdf_adj_init = betta*c1(T_sim)^(1-eps).*cprime_init.^(eps-1)*( 1 / E_reverse_init )^(1-gama-eps).*cprime_init.^(1-gama-eps).*pi_prime.^(-1);

        Resid_init(k_init) = abs(1 - R_guess*prob*sdf_adj_init);

        if k_init > 1 && Resid_init(k_init) > Resid_init(k_init-1)
            break 
        end 

    end

    [Resid_init,Index_init] = min(Resid_init);
    R_ann = R_guess_stack_init(Index_init);
    R(T_sim) = R_ann;

    %-------------
    %Initial old
    %------------- 

    cprime_init = alfa*(1+n)*y_prime + R_ann*bstar.*pi_prime.^(-1);
    c_power_init = cprime_init.^(1-gama);
    E_c_power_init = prob*c_power_init;
    Utility_init =  1/(1-gama)*( c1(T_sim)^eps + betta*( E_c_power_init )^(eps/(1-gama)) )^((1-gama)/eps);

    else

    %-------------
    %Initial old
    %------------- 
    pi_prime_init = pistar*exp(e_vec(T_sim+1))/(y(T_sim+1)/y(T_sim))^(1-dummy_IT);
    cprime_init = alfa*(1+n)*y(T_sim+1) + R(T_sim)*bstar/pi_prime_init;
    Utility_init = 1/(1-gama)*( c1(T_sim)^eps + betta*( cprime_init^(1-gama) )^(eps/(1-gama)) )^((1-gama)/eps);

    %Utility_init =  Utility(T_sim);  %if use expected U of initial old

    end

    if dummy_IT == 1
        U_IT_init = Utility_init;
    else
        load("test.mat","U_IT","U_IT_init")
        %Utility_init = U_IT_init;
    end

%------------
%Reform
%------------

for t=T_sim+1:T_fin

    R_lag = R(t-1);
    y_lag = y(t-1);

    pi(t) = pistar*exp(e_vec(t))/(y(t)/y_lag)^(1-dummy_IT);
    tau(t) = ( gbar + (R_lag/((1+n)*pi(t)) - 1)*bstar ) / y(t);
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

    EV(t) = prob*c_power;

    Utility(t) =  1/(1-gama)*( c1(t)^eps + betta*( EV(t) )^(eps/(1-gama)) )^((1-gama)/eps);

    if t==T_sim+1
        U_sum(t) = omega^(t-1-T_sim)*Utility(t) + omega^(-1)*Utility_init;
    else
        U_sum(t) = omega^(t-1-T_sim)*Utility(t);
    end

    if dummy_IT == 0
        Lambda_NIT(t) = 100*( (Utility(t)/U_IT(t))^(1/(1-gama)) - 1);  %Consumption equiv. welfare gain
        Lambda_init = 100*( (Utility_init/U_IT_init)^(1/(1-gama)) - 1);
    end

end

U_sum_tot = sum(U_sum(T_sim+1:T_fin));

if dummy_IT == 1
    U_IT = Utility;
    U_sum_IT = U_sum_tot;
end





