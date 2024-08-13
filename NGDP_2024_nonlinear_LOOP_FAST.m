%NGDP_2024_nonlinear_LOOP_FAST
%Nonlinear simulation of the simple model (Algorithm in Supp Appendix),
%this faster version of the code avoids stacking over many simulations. 
%Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). Any errors are my own.

%clear

alfa = 0.3;  
betta = 0.85;  
chi = 1;
dummy_IT = 1;
gama = 5;  
eps = 0.5;  
n = 0.4;  
pistar = 1.8;
ybar = 1;
gbar = 0.15;  
phi = 0.5;
sig_e = 0.025;  
sig_A = 0.05; 

n_sim = 100;
T_sim = 310; drop = 10;
T_sim2 = T_sim - drop;

N_guess0 = 610;
N_guess = 800;
y_init = 1;

n_loop = 15;
chi_stack = linspace(0.9,1.045,n_loop);
R_lower = -1.05; R_upper = 1.05;
R_l = -0.085; R_u = 0.085;

n_states = 5;  %No. of states
prob = ones(1,n_states^2); prob = prob / sum(prob);

%Shocks
sigma = sig_e;
Discretization_short
y1 = e_i;
sigma = sig_A;
Discretization_short
x1 = e_i;
%Grid of values
[X,Y] = meshgrid(x1,y1);  states = [Y(:) X(:)];
y_prime = ybar*exp(states(:,2));

Stack_utility = NaN(n_loop,1); Stack_tau = Stack_utility; Stack_ftau = Stack_utility; Stack_b = Stack_utility; 
Stack_c1 = Stack_utility; Stack_c2 = Stack_utility; var_tau = Stack_b; var_ftau = Stack_b; Lambda = Stack_utility;
U_IT = Stack_b; tau_IT = U_IT; ftau_IT = U_IT; vtau_IT = U_IT; c1_IT = U_IT; c2_IT = U_IT; vc1_IT = U_IT; vc2_IT = U_IT; 
var_c1 = Stack_utility; var_c2 = Stack_utility; U_mean_NIT = var_c1; ftau_mean = var_c1;

Resid_max = U_IT; Lambda_NIT = Resid_max; Index_mini = Resid_max; Index_maxi = Resid_max; Resid_max2 = Resid_max;  
Lambda_tau = Resid_max; Lambda_mean_NIT = Resid_max; Lambda_var_NIT = Resid_max; U_mean_IT = Resid_max;

for m=1:n_loop

    chi = chi_stack(m);
    NGDP_steady_state

    bstar = b_root;
    Rstar = pistar*chi*(1+n);
    Stack_b(m) = bstar;

    R_init = Rstar; 

    %Guesses for interest rate
    R_guess_stack0 = R_init + linspace(R_lower,R_upper,N_guess0);
    %R_guess_stack = R_guess_stack(randperm(N_guess));

    %---------------------------
    %Stochastic simulations
    %---------------------------
    Dum = zeros(n_sim,1); Max_resid = Dum; Index_loc = Dum; Index_loc2 = Dum;  
    c1_vec = Dum; c2_vec =  Dum; tau_vec =  Dum; ftau_vec =  Dum; Utility_vec =  Dum;
    c1_vec2 = Dum; c2_vec2 =  Dum; tau_vec2 =  Dum; ftau_vec2 =  Dum; 

    for j=1:n_sim

        rng(500+j)
        y = ybar*exp(randn(T_sim,1)*sig_A);
        e_vec = randn(T_sim,1)*sig_e;

        NGDP_2024_nonlinear_SIM

        Max_resid(j) = max(Max_Resid);
        Index_loc(j) = max(Index);
        Index_loc2(j) = min(Index);
    
        %Store and stack key variables

        %For means
        c1_vec(j) = sum(c1(drop+1:end));
        c2_vec(j) = sum(c2(drop+1:end));
        tau_vec(j) = sum(tau(drop+1:end));
        ftau_vec(j) = sum(ftau(drop+1:end));
        Utility_vec(j) = sum(Utility(drop+1:end));

        %For variances
        c1_vec2(j) = sum(c1(drop+1:end).^2);
        c2_vec2(j) = sum(c2(drop+1:end).^2);
        tau_vec2(j) = sum(tau(drop+1:end).^2);
        ftau_vec2(j) = sum(ftau(drop+1:end).^2);

    end

    Resid_max(m) = max(Max_resid);
    Resid_max2(m) = max(Resid_check);

    Index_mini(m) = min(Index_loc2);
    Index_maxi(m) = max(Index_loc);

    Stack_utility(m) = sum(Utility_vec)/(n_sim*T_sim2);
    Stack_c1(m) = sum(c1_vec)/(n_sim*T_sim2);
    Stack_c2(m) = sum(c2_vec)/(n_sim*T_sim2);
    Stack_tau(m) = sum(tau_vec)/(n_sim*T_sim2);
    Stack_ftau(m) = sum(ftau_vec)/(n_sim*T_sim2);
    
    var_c1(m) = sum(c1_vec2)/(n_sim*T_sim2) - (Stack_c1(m))^2;
    var_c2(m) = sum(c2_vec2)/(n_sim*T_sim2) - (Stack_c2(m))^2;
    var_tau(m) = sum(tau_vec2)/(n_sim*T_sim2) - (Stack_tau(m))^2;
    var_ftau(m) = sum(ftau_vec2)/(n_sim*T_sim2) - (Stack_ftau(m))^2;

end

    if dummy_IT == 1
        U_IT = Stack_utility;
        c1_IT = Stack_c1;
        c2_IT = Stack_c2;
        vc1_IT = var_c1;
        vc2_IT = var_c2;
        tau_IT = Stack_tau;
        ftau_IT = Stack_ftau;
        vtau_IT = var_tau;

        Lambda = 100*( (Stack_utility.*1/U_IT(1)).^(1/(1-gama)) - 1);  %Consumption equiv. welfare gain
        for i=1:n_loop
            U_mean_IT(i) = (1/(1-gama))*( c1_IT(i)^eps + betta*c2_IT(i)^eps  )^((1-gama)/eps);
        end

        save sim_IT.mat

        figure(1)
        hold on, subplot(2,2,3),  plot(Stack_b, Lambda, 'k','DisplayName', 'IT', 'LineWidth', 1), title('Welfare gain vs bond supply: Stochastic case'), hold on,
        xlabel('Steady state bond supply'), ylabel('% c.e. welfare gain')

        h = gcf;  %gcf - get currrent figure - save the figure handle for further modification
        savefig(h, 'first_run_plot_2.fig')

        save sim_nonlin_IT.mat

    elseif dummy_IT == 0

        load("sim_nonlin_IT.mat","U_IT","c1_IT","c2_IT","vc1_IT","vc2_IT","tau_IT","ftau_IT","U_mean_IT")
        
        Lambda = 100*( (Stack_utility.*1/U_IT(1)).^(1/(1-gama)) - 1);  %Consumption equiv. welfare gain
        Lambda_NIT = 100*( (Stack_utility.*U_IT.^(-1)).^(1/(1-gama)) - 1);  %Consumption equiv. welfare gain
        
        c1_NIT = Stack_c1.*c1_IT.^(-1); c2_NIT = Stack_c2.*c2_IT.^(-1); 
        vc1_NIT = var_c1.*vc1_IT.^(-1); vc2_NIT = var_c2.*vc2_IT.^(-1);
        
        tau_NIT = Stack_tau.*tau_IT.^(-1); vtau_NIT = var_tau.*vtau_IT.^(-1);
 
        for i=1:n_loop
        
            Lambda_tau(i) = 100*(sqrt(Stack_ftau(i)/ftau_IT(i)) -1);
            ftau_mean(i) = 100*(sqrt(Stack_tau(i)^2/tau_IT(i)^2) -1);
            U_mean_NIT(i) = (1/(1-gama))*( Stack_c1(i)^eps + betta*Stack_c2(i)^eps  )^((1-gama)/eps);
            Lambda_mean_NIT(i) = 100*( (U_mean_NIT(i)/U_mean_IT(i))^(1/(1-gama)) - 1);
            Lambda_var_NIT(i) = Lambda_NIT(i) - Lambda_mean_NIT(i);

        end

        ftau_var = Lambda_tau - ftau_mean;

          
        h = openfig('first_run_plot_2.fig', 'reuse');    % open the previously saved figure and add new plots
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

Index_lower = min(Index_mini)
Index_upper = max(Index_maxi)
max(Resid_max)
max(Resid_max2)




