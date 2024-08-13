%NGDP_KL2_Plotter
%Written by M. Hatcher based on the code by M. Hatcher and J. Lyu

%--------------
%Plot results
%--------------
Lambda = NaN(n_loop,1);

if dummy_IT==1
    U_mean_IT =  NaN(n_loop,1);

    U_IT = Stack_utility;  
    c1_IT = Stack_c1; c2_IT = Stack_c2; vc1_IT = var_c1; vc2_IT = var_c2; comp_IT = Stack_comp;
    tau_IT = Stack_tau; vtau_IT = var_tau; r_IT = Stack_r; vr_IT = var_r; y_IT = Stack_y; 
    k_IT = Stack_k; l_IT = Stack_l; rk_IT = Stack_rk; vrk_IT = var_rk;

    for i=1:n_loop
        Lambda(i) = 100*( (Stack_utility(i)/Stack_utility(1))^(1/thetta*(1-gama)) - 1);  %Consumption equiv. welfare gain 
        U_mean_IT(i) = (1/(1-gama))*( comp_IT(i)^eps + betta*c2_IT(i)^(thetta*eps)  )^((1-gama)/eps);
    end

    save sim_IT_KL2.mat

else

Lambda_NIT = Lambda; c1_NIT = Lambda; c2_NIT = Lambda; vc1_NIT = Lambda; vc2_NIT = Lambda; tau_NIT = Lambda; vtau_NIT = Lambda;
r_NIT = Lambda; vr_NIT = Lambda; y_NIT = Lambda; k_NIT = Lambda; l_NIT = Lambda; rk_NIT = Lambda; vrk_NIT = Lambda;
U_mean_NIT = Lambda; Lambda_mean_NIT = Lambda; Lambda_var_NIT = Lambda;

load("sim_IT_KL2.mat","U_IT","c1_IT","c2_IT","vc1_IT","vc2_IT","tau_IT","vtau_IT","r_IT","vr_IT","y_IT","k_IT","l_IT","rk_IT","vrk_IT")

    for i=1:n_loop
        Lambda(i) = 100*( (Stack_utility(i)/U_IT(1))^(1/thetta*(1-gama)) - 1);  %Consumption equiv. welfare gain
        Lambda_NIT(i) = 100*( (Stack_utility(i)/U_IT(i))^(1/thetta*(1-gama)) - 1);  %Consumption equiv. welfare gain
        c1_NIT(i) = Stack_c1(i)/c1_IT(i); c2_NIT(i) = Stack_c2(i)/c2_IT(i); 
        vc1_NIT(i) = var_c1(i)/vc1_IT(i); vc2_NIT(i) = var_c2(i)/vc2_IT(i);
        
        tau_NIT(i) = Stack_tau(i)/tau_IT(i); vtau_NIT(i) = var_tau(i)/vtau_IT(i);
        r_NIT(i) = Stack_r(i)/r_IT(i);
        vr_NIT(i) = var_r(i)/vr_IT(i);

        y_NIT(i) = Stack_y(i)/y_IT(i); k_NIT(i) = Stack_k(i)/k_IT(i); l_NIT(i) = Stack_l(i)/l_IT(i); 
        rk_NIT(i) = Stack_rk(i)/rk_IT(i); vrk_NIT(i) = var_rk(i)/vrk_IT(i);
        
        U_mean_NIT(i) = (1/(1-gama))*( Stack_comp(i)^eps + betta*Stack_c2(i)^(thetta*eps)  )^((1-gama)/eps);
        Lambda_mean_NIT(i) = 100*( (U_mean_NIT(i)/U_mean_IT(i))^(1/thetta*(1-gama)) - 1);
        Lambda_var_NIT(i) = Lambda_NIT(i) - Lambda_mean_NIT(i);

    end

    figure(1)
    subplot(2,3,1), plot(Stack_b, c2_NIT, '--k', 'LineWidth', 1), hold on, plot(Stack_b, c1_NIT, 'k', 'LineWidth', 1) 
    title('Consumption means'), hold on, xlabel('Fixed bond supply'), ylabel('Ratio')
    subplot(2,3,2), plot(Stack_b, vc2_NIT, '--k', 'LineWidth', 1), hold on, plot(Stack_b, vc1_NIT, 'k', 'LineWidth', 1) 
    title('Consumption variances'), hold on, xlabel('Fixed bond supply'), ylabel('Ratio')
    subplot(2,3,3), plot(Stack_b, y_NIT, 'k', 'LineWidth', 1),  
    title('Mean output: $E[y]$'), hold on, xlabel('Fixed bond supply'), ylabel('Ratio')
    subplot(2,3,4), plot(Stack_b, k_NIT, 'k', 'LineWidth', 1) 
    title('Capital: $E[k]$'), hold on, xlabel('Fixed bond supply'), ylabel('Ratio')
    subplot(2,3,5), plot(Stack_b, r_NIT, 'k', 'LineWidth', 1), hold on, plot(Stack_b, rk_NIT, '--k', 'LineWidth', 1)
    title('Expected returns: $E[r_t^i]$'), hold on, xlabel('Fixed bond supply'), ylabel('Ratio')
    subplot(2,3,6), plot(Stack_b, vr_NIT, 'k', 'LineWidth', 1), hold on, plot(Stack_b, vrk_NIT, '--k', 'LineWidth', 1)  
    title('Return risk: $var[r_t^i]$'), hold on, xlabel('Fixed bond supply'), ylabel('Ratio')

    figure(2)
    subplot(1,2,1), plot(Stack_b,Stack_Lambda_ss,'k','LineWidth', 1), hold on, plot(Stack_b,Stack_Lambda_ss1,'--k','LineWidth', 1)
    title('Welfare'), hold on, xlabel('Steady state bond supply'), ylabel('% cons. equiv.'), hold on,
    subplot(1,2,2), plot(Stack_b,Tax_burden_ss,'k','LineWidth', 1)
    title('Tax burden'), hold on, xlabel('Steady state bond supply'), ylabel('% cons. equiv.')

end

NGDP_KL2_Run_1







