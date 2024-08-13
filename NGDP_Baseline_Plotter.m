%NGDP_Baseline_Plotter
%Written by M. Hatcher and J. Lyu
%Used in NGDP_2024_Baseline.mod and its variants

%--------------
%Plot results
%--------------
Lambda = NaN(n_loop,1); 

if dummy_IT==1

    U_mean_IT = NaN(n_loop,1);

    U_IT = Stack_utility;   %Record and store IT results
    c1_IT = Stack_c1; c2_IT = Stack_c2; vc1_IT = var_c1; vc2_IT = var_c2;
    tau_IT = Stack_tau; vtau_IT = var_tau; ftau_IT = Stack_ftau;
    r_IT = Stack_r; vr_IT = var_r;

    for i=1:n_loop
        Lambda(i) = 100*( (Stack_utility(i)/Stack_utility(1))^(1/(1-gama)) - 1);  %Consumption equiv. welfare gain 
        U_mean_IT(i) = (1/(1-gama))*( c1_IT(i)^eps + betta*c2_IT(i)^eps  )^((1-gama)/eps);
    end

    save sim_IT.mat

else

load("sim_IT.mat","U_IT","c1_IT","c2_IT","vc1_IT","vc2_IT","tau_IT","vtau_IT","r_IT","vr_IT")

Lambda_NIT = NaN(n_loop,1); c1_NIT = Lambda_NIT; c2_NIT = c1_NIT; vc1_NIT = c1_NIT; vc2_NIT = c1_NIT; tau_NIT = c1_NIT; 
vtau_NIT = c1_NIT; Lambda_tau = c1_NIT; ftau_mean = c1_NIT; ftau_var = c1_NIT; r_NIT = c1_NIT; vr_NIT = c1_NIT;
U_mean_NIT = c1_NIT; Lambda_mean_NIT = c1_NIT; Lambda_var_NIT = c1_NIT;

    for i=1:n_loop
        Lambda(i) = 100*( (Stack_utility(i)/U_IT(1))^(1/(1-gama)) - 1);  %Consumption equiv. welfare gain
        Lambda_NIT(i) = 100*( (Stack_utility(i)/U_IT(i))^(1/(1-gama)) - 1);  %Consumption equiv. welfare gain
        
        c1_NIT(i) = Stack_c1(i)/c1_IT(i); c2_NIT(i) = Stack_c2(i)/c2_IT(i); 
        vc1_NIT(i) = var_c1(i)/vc1_IT(i); vc2_NIT(i) = var_c2(i)/vc2_IT(i);
        
        tau_NIT(i) = Stack_tau(i)/tau_IT(i); vtau_NIT(i) = var_tau(i)/vtau_IT(i);
        Lambda_tau(i) = 100*(sqrt(Stack_ftau(i)/ftau_IT(i)) -1);
        ftau_mean(i) = 100*(sqrt(Stack_tau(i)^2/tau_IT(i)^2) -1);
        ftau_var(i) = Lambda_tau(i) - ftau_mean(i);
        r_NIT(i) = Stack_r(i)/r_IT(i);
        vr_NIT(i) = var_r(i)/vr_IT(i);
        
        U_mean_NIT(i) = (1/(1-gama))*( Stack_c1(i)^eps + betta*Stack_c2(i)^eps  )^((1-gama)/eps);
        Lambda_mean_NIT(i) = 100*( (U_mean_NIT(i)/U_mean_IT(i))^(1/(1-gama)) - 1);
        Lambda_var_NIT(i) = Lambda_NIT(i) - Lambda_mean_NIT(i);

    end

    figure(1)
    subplot(2,3,1), plot(Stack_b, c1_NIT, 'k', 'LineWidth', 1), hold on, plot(Stack_b, c2_NIT, '--k', 'LineWidth', 1) 
    title('Consumption means'), hold on, xlabel('Fixed bond supply'), ylabel('Ratio')
    subplot(2,3,2), plot(Stack_b, vc1_NIT, 'k', 'LineWidth', 1), hold on, plot(Stack_b, vc2_NIT, '--k', 'LineWidth', 1) 
    title('Consumption variances'), hold on, xlabel('Fixed bond supply'), ylabel('Ratio')
    subplot(2,3,3), plot(Stack_b, tau_NIT, 'k', 'LineWidth', 1),  
    title('Mean taxes: $E[\tau_t]$'), hold on, xlabel('Fixed bond supply'), ylabel('Ratio')
    subplot(2,3,4), plot(Stack_b, vtau_NIT, 'k', 'LineWidth', 1), 
    title('Variance of taxes: $var[\tau_t]$'), hold on, xlabel('Fixed bond supply'), ylabel('Ratio')
    subplot(2,3,5), plot(Stack_b, r_NIT, 'k', 'LineWidth', 1),
    title('Expected return: E[r_t^n]$'), hold on, xlabel('Fixed bond supply'), ylabel('Ratio')
    subplot(2,3,6), plot(Stack_b, vr_NIT, 'k', 'LineWidth', 1),  
    title('Return risk: $var[r_t^n]$'), hold on, xlabel('Fixed bond supply'), ylabel('Ratio')

end

if phi == 0.5
    %------------------Generate the initial plot (1st run) ------------
    NGDP_Baseline_Run_1
    U_init = U_ss(1);
    save sim_phi_0.5.mat

elseif phi == 0 
    %-------------------- Adding the plot (lump-sum taxes) ----------
    load("sim_phi_0.5.mat","U_init")
    NGDP_Baseline_Run_2
end







