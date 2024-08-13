%NGDP_2024_nonlinear_transition_LOOP
%Nonlinear simulation of transition to NIT in the simple model (Algorithm
%in Supp Appendix) for realized shocks, average utility over many simumations. 
%Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). Any errors are my own.

%clear

alfa = 0.3;  
betta = 0.85;  
dummy_IT = 1;
dummy_NIT = 1 - dummy_IT;
gama = 5;  
eps = 0.5;  
n = 0.4;  
pistar = 1.8;
ybar = 1;
gbar = 0.15;  
phi = 0.5;
sig_e = 0.025;  
sig_A = 0.05; 

omega = 0.95;  %for social welfare analysis
T_sim = 5;  %5  
T_fin = 255;  %10  %255
n_sim = 1100; 

N_guess0 = 800;
N_guess = 2000;
N_guess_init = 2000;
y_init = 1;

Announced = 1;
R_lower = -1.1; R_upper = 1.1;
R_l = -0.1; R_u = 0.1;

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

chi = 0.9659;
NGDP_steady_state

bstar = b_root;  Rstar = pistar*chi*(1+n);

R_init = Rstar; 
R_guess_stack0 = R_init + linspace(R_lower,R_upper,N_guess0);
R_guess_stack_init = R_init + linspace(R_lower,R_upper,N_guess_init);

%---------------------------
%Stochastic simulations
%--------------------------- 
T_sim2 = T_fin - T_sim; Index_loc = NaN(n_sim,1);  Index_loc2 = Index_loc;
Max_resid = NaN(n_sim,1); Resid_init_stack = Max_resid; U_init_stack = NaN(n_sim,1); 
U_stack = NaN(n_sim,T_sim2); U_sum_stack = U_stack; U_sum_stack0 = U_init_stack; Resid_init = NaN;

for j = 1:n_sim

    rng(5E5+j)
    y = ybar*exp(randn(T_fin,1)*sig_A);
    e_vec = randn(T_fin,1)*sig_e;

    NGDP_2024_nonlinear_SIM_transition
        
    Index_loc(j) = max(Index);
    Index_loc2(j) = min(Index);
    Max_resid(j) = max(Max_Resid);

    Resid_init_stack(j) = Resid_init;

    U_stack(j,:) = Utility(T_sim+1:T_fin)';
    U_init_stack(j) = Utility_init;
    U_sum_stack(j,:) = U_sum(T_sim+1:T_fin)';
    U_sum_stack0(j) = sum(U_sum_stack(j,:));

end  

if dummy_IT == 1
        Ue_IT = mean(U_stack); Ue_IT_init = mean(U_init_stack);
        Ue_IT = [Ue_IT_init Ue_IT]';
        SW_IT = (1-omega)*mean(U_sum_stack0);
        save test.mat
else
        load("test.mat","U_IT","U_IT_init","SW_IT")
        Ue_NIT = mean(U_stack); Ue_NIT_init = mean(U_init_stack);
        Ue_NIT = [Ue_NIT_init Ue_NIT]';
        SW = (1-omega)*mean(U_sum_stack0);
        Lambda_SW = 100*( (SW/SW_IT)^(1/(1-gama)) - 1)  %Consumption equiv. social welfare gain

    Lambda_NIT_plot = 100*( (Ue_NIT.*Ue_IT.^(-1)).^(1/(1-gama)) - 1);  %Consumption equiv. welfare gain;
    Periods = T_sim:T_fin; Periods = Periods - T_sim;
    hold on, plot(Periods,Lambda_NIT_plot,'--k','LineWidth', 1),
    xlabel('Generations (d.o.b.)'), ylabel('% c.e. welfare gain')
end

max(Max_resid)

if Announced == 1
    Resid_init_max = max(Resid_init_stack)
end

