%NGDP_Baseline_Run_2

h = openfig('first_run_plot.fig', 'reuse');    % open the previously saved figure and add new plots
figure(h)  % Make sure it's the current figure
Lambda_ss = NaN(n_loop,1);

for i=1:n_loop
    Lambda_ss(i) = 100*( (U_ss(i)/U_init)^(1/(1-gama))  - 1);   %Consumption equiv. welfare gain 
end

hold on, subplot(1,2,1), plot(r_ss, b_ss, '--k','DisplayName', '\Phi = 0', 'LineWidth', 1), hold on, 
title('Bond supply vs interest rate'), xlabel('Steady state interest rate'), ylabel('Steady state bond supply')
hold on, subplot(1,2,2), xline( 1.4,'--k' ), hold on, plot(r_ss, Lambda_ss, '--k','LineWidth', 1) 
title('Interest rate vs steady state welfare'), xlabel('Steady state interest rate'), ylabel('% c.e. welfare gain'); % Add data from the second Dynare run
legend('show');  %Refresh the legend to include the new line
savefig(h, 'combined_run_plot.fig');

if dummy_IT==0

    U_NIT = Stack_utility;

    load("sim_IT.mat","U_IT")
    Lambda_NIT = NaN(n_loop,1); Lambda_var_NIT = Lambda_NIT;

    for i=1:n_loop
        Lambda_NIT(i) = 100*( (U_NIT(i)/U_IT(i))^(1/(1-gama)) - 1);  %Consumption equiv. welfare gain 
        Lambda_var_NIT(i) = Lambda_NIT(i) - Lambda_mean_NIT(i);
    end

    h = openfig('first_run_plot_1.fig', 'reuse');    % open the previously saved figure and add new plots
    figure(h)  % Make sure it's the current figure

    subplot(2,2,1), hold on, plot(Stack_b, Lambda_tau, 'k', 'LineWidth', 1), 
    title('Tax burden: increase under NIT'), hold on, xlabel('Fixed bond supply'), ylabel('% tax equiv.')
    subplot(2,2,2), plot(Stack_b, ftau_mean, 'k', 'LineWidth', 1), hold on, plot(Stack_b, ftau_var, '--k', 'LineWidth', 1)
    title('Tax burden: decomposition'), hold on, xlabel('Fixed bond supply'), ylabel('% tax equiv.')
    subplot(2,2,3), plot(Stack_b, Lambda, '--k','DisplayName', 'IT', 'LineWidth', 1), title('Welfare gain vs bond supply: Stochastic case'), hold on,
    title('Welfare'), hold on, xlabel('Fixed bond supply'), ylabel('% cons. equiv.')
    subplot(2,2,4), plot(Stack_b, Lambda_NIT, 'k', 'LineWidth', 2), hold on, plot(Stack_b, Lambda_mean_NIT, 'k', 'LineWidth', 1), hold on, plot(Stack_b, Lambda_var_NIT, '--k', 'LineWidth', 1)
    title('Welfare loss decomposition'), hold on, xlabel('Fixed bond supply'), ylabel('% cons. equiv.')

end