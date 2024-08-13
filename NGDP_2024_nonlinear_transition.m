%NGDP_2024_nonlinear_transition
%Nonlinear simulation of transition to NIT in the simple model (Algorithm
%in Supp Appendix) for zero realized shocks, inc. social welfare computation. 
%Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). Any errors are my own.

%clear

alfa = 0.3;  
betta = 0.85;  
dummy_IT = 0;
dummy_NIT = 1-dummy_IT;
gama = 5;  
eps = 0.5;  
n = 0.4;  
pistar = 1.8;
ybar = 1;
gbar = 0.15;  
phi = 0.5;
sig_e = 0.025;  
sig_A = 0.05; 

omega = 0.5;  %for welfare analysis
T_sim = 5;    %5 
T_fin = 10;   %255 %10

N_guess0 = 800;
N_guess = 2000;
N_guess_init = 2000;
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
taustar = ( gbar + (chi-1)*bstar ) / ybar;
c1star = (1-alfa-taustar)*ybar - phi*taustar^2 - bstar;
c2star = alfa*(1+n)*ybar + chi*(1+n)*bstar;
utilitystar = (1/(1-gama))*(c1star^eps + betta*c2star^eps)^((1-gama)/eps);

R_init = Rstar; 
R_guess_stack0 = R_init + linspace(R_lower,R_upper,N_guess0);
R_guess_stack_init = R_init + linspace(R_lower,R_upper,N_guess);

%---------------------------
%Stochastic simulations
%---------------------------    
rng(1)
y = ybar*exp(randn(T_fin,1)*sig_A*0);
e_vec = randn(T_fin,1)*sig_e*0;

NGDP_2024_nonlinear_SIM_transition

if dummy_IT == 1
    SW_IT = (1-omega)*U_sum_IT;
    save test.mat
else
    SW = (1-omega)*U_sum_tot;
    Lambda_SW = 100*( (SW/SW_IT)^(1/(1-gama)) - 1)  %Consumption equiv. social welfare gain
    Lambda_NIT_plot = [Lambda_init; Lambda_NIT(T_sim+1:T_fin)];
    Periods = T_sim:T_fin; Periods = Periods - T_sim;
    hold on, plot(Periods,Lambda_NIT_plot,'--k','LineWidth', 1),
    xlabel('Generations (d.o.b.)'), ylabel('% c.e. welfare gain'), %title('Unannounced')
end

Index_loc = max(Index)
Index_loc2 = min(Index)
Max_resid = max(Max_Resid)

