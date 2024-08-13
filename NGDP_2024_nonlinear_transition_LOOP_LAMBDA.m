%NGDP_20204_nonlinear_transition

%clear

alfa = 0.3;  
betta = 0.85;  
dummy_IT = 0;
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

%Social discount factors
omega = 0.975;  %for social welfare analysis
omega_1 = 0.95;  %for social welfare analysis
omega_2 = 0.90;  %for social welfare analysis
omega_3 = 0.75;  %for social welfare analysis
omega_4 = 0.50;  %for social welfare analysis

T_sim = 5;  %5  %500
T_fin = 255;  %10  %255
n_sim = 1100;  

N_guess0 = 800;
N_guess = 800;  %for speed in social welfare analysis
N_guess_init = 800;   %for speed in social welfare analysis
y_init = 1;

Announced = 0;
R_lower = -1.1; R_upper = 1.1;
R_l = -0.1; R_u = 0.1;

n_states = 5;  %No. of states
prob = ones(1,n_states^2); prob = prob / sum(prob);
Resid_init = NaN(N_guess,1);

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
Max_resid = NaN(n_sim,1); Resid_init_stack = Max_resid; U_init_stack = NaN(n_sim,1); Max_resid2 = Max_resid;
U_stack = NaN(n_sim,T_sim2); U_sum_stack = U_stack; U_sum_stack_1 = U_stack; U_sum_stack_2 = U_stack;
U_sum_stack_3 = U_stack; U_sum_stack_4 = U_stack; U_sum_stack0 = U_init_stack; U_sum_stack1 = U_init_stack;
U_sum_stack2 = U_init_stack;U_sum_stack3 = U_init_stack; U_sum_stack4 = U_init_stack; Resid_init = NaN;

for j = 1:n_sim

    rng(5E5+j)
    y = ybar*exp(randn(T_fin,1)*sig_A);
    e_vec = randn(T_fin,1)*sig_e;

    NGDP_2024_nonlinear_SIM_transition_LAMBDA
        
    Index_loc(j) = max(Index);
    Index_loc2(j) = min(Index);

    Max_resid(j) = max(Max_Resid);
    Max_resid2(j) = max(Resid_check);

    Resid_init_stack(j) = Resid_init;

    U_stack(j,:) = Utility(T_sim+1:T_fin)';
    U_init_stack(j) = Utility_init;

    U_sum_stack0(j) = U_sum_tot;
    U_sum_stack1(j) = U_sum_tot_1;
    U_sum_stack2(j) = U_sum_tot_2;
    U_sum_stack3(j) = U_sum_tot_3;
    U_sum_stack4(j) = U_sum_tot_4;

end  

if dummy_IT == 1
        Ue_IT = mean(U_stack); Ue_IT_init = mean(U_init_stack);
        Ue_IT = [Ue_IT_init Ue_IT]';

        SW_IT = (1-omega)*mean(U_sum_stack0);
        SW_IT_1 = (1-omega_1)*mean(U_sum_stack1);
        SW_IT_2 = (1-omega_2)*mean(U_sum_stack2);
        SW_IT_3 = (1-omega_3)*mean(U_sum_stack3);
        SW_IT_4 = (1-omega_4)*mean(U_sum_stack4);
        
        save test_lambda.mat
else
        load("test_lambda.mat","U_IT","U_IT_init","SW_IT","SW_IT_1","SW_IT_2","SW_IT_3","SW_IT_4")
        Ue_NIT = mean(U_stack); Ue_NIT_init = mean(U_init_stack);
        Ue_NIT = [Ue_NIT_init Ue_NIT]';

        SW = (1-omega)*mean(U_sum_stack0);
        SW_1 = (1-omega_1)*mean(U_sum_stack1);
        SW_2 = (1-omega_2)*mean(U_sum_stack2);
        SW_3 = (1-omega_3)*mean(U_sum_stack3);
        SW_4 = (1-omega_4)*mean(U_sum_stack4);

        Lambda_SW = 100*( (SW/SW_IT)^(1/(1-gama)) - 1)  %Consumption equiv. social welfare gain
        Lambda_SW_1 = 100*( (SW_1/SW_IT_1)^(1/(1-gama)) - 1)  %Consumption equiv. social welfare gain
        Lambda_SW_2 = 100*( (SW_2/SW_IT_2)^(1/(1-gama)) - 1)  %Consumption equiv. social welfare gain
        Lambda_SW_3 = 100*( (SW_3/SW_IT_3)^(1/(1-gama)) - 1)  %Consumption equiv. social welfare gain
        Lambda_SW_4 = 100*( (SW_4/SW_IT_4)^(1/(1-gama)) - 1)  %Consumption equiv. social welfare gain

    Lambda_NIT_plot = 100*( (Ue_NIT.*Ue_IT.^(-1)).^(1/(1-gama)) - 1);  %Consumption equiv. welfare gain;
    Periods = T_sim:T_fin; Periods = Periods - T_sim;
    hold on, plot(Periods,Lambda_NIT_plot,'--k','LineWidth', 1),
    xlabel('Generations (d.o.b.)'), ylabel('% c.e. welfare gain')
end

max(Max_resid)
max(Max_resid2)

if Announced == 1
    Resid_init_max = max(Resid_init_stack)
end

