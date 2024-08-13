%NGDP_Baseline_Run_1

figure(3)
hold on, subplot(1,2,1), plot(r_ss, b_ss, 'k','DisplayName', '\Phi = 0.5', 'LineWidth', 1), hold on, 
title('Bond supply vs interest rate'), xlabel('Steady state interest rate'), ylabel('Steady state bond supply')
hold on, subplot(1,2,2), plot(r_ss, Lambda_ss, 'k','LineWidth', 1) 
title('Interest rate vs steady state welfare'), xlabel('Steady state interest rate'), ylabel('% c.e. welfare gain')
%hold on, subplot(1,3,3), plot(b_ss, Lambda_ss, 'k') 
%title('Interest rate versus steady state welfare'), xlabel('Steady state interest rate'), ylabel('% c.e. welfare gain')

h = gcf;  %gcf - get currrent figure - save the figure handle for further modification
savefig(h, 'first_run_plot.fig')

if dummy_IT==1
    figure(4)
    hold on, subplot(2,2,3),  plot(Stack_b, Lambda, 'k','DisplayName', 'IT', 'LineWidth', 1), title('Welfare gain vs bond supply: Stochastic case'), hold on,
    xlabel('Steady state bond supply'), ylabel('% c.e. welfare gain')

    U_IT = Stack_utility;
    save sim_IT.mat

    h = gcf;  %gcf - get currrent figure - save the figure handle for further modification
    savefig(h, 'first_run_plot_1.fig')

elseif dummy_IT==0

    U_NIT = Stack_utility;
    load("sim_IT.mat","U_IT")
    Lambda_NIT = NaN(n_loop,1);

    for i=1:n_loop
        Lambda_NIT(i) = 100*( (U_NIT(i)/U_IT(i))^(1/(1-gama)) - 1);  %Consumption equiv. welfare gain 
    end

    h = openfig('first_run_plot_1.fig', 'reuse');    % open the previously saved figure and add new plots
    figure(h)  % Make sure it's the current figure

    subplot(2,2,1), hold on, plot(Stack_b, Lambda_tau, 'k', 'LineWidth', 1), 
    title('Tax burden: increase under NIT'), hold on, xlabel('Fixed bond supply'), ylabel('% tax equiv.')
    subplot(2,2,2), plot(Stack_b, ftau_mean, 'k', 'LineWidth', 1), hold on, plot(Stack_b, ftau_var, '--k', 'LineWidth', 1)
    title('Tax burden: decomposition'), hold on, xlabel('Fixed bond supply'), ylabel('% tax equiv.')
    subplot(2,2,3), plot(Stack_b, Lambda, '--k','DisplayName', 'NIT', 'LineWidth', 1), title('Welfare gain vs bond supply: Stochastic case'), hold on,
    title('Welfare'), hold on, xlabel('Fixed bond supply'), ylabel('% cons. equiv.')
    subplot(2,2,4), plot(Stack_b, Lambda_NIT, 'k', 'LineWidth', 2), hold on, plot(Stack_b, Lambda_mean_NIT, 'k', 'LineWidth', 1), hold on, plot(Stack_b, Lambda_var_NIT, '--k', 'LineWidth', 1)
    title('Welfare loss decomposition'), hold on, xlabel('Fixed bond supply'), ylabel('% cons. equiv.')
 

end