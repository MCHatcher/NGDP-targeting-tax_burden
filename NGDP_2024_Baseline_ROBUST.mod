//Code for robsutness checks in the simple model (Supplementary Appendix). 
//Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). Any errors are my own.

//-----------------------------------------
//1. Variable declaration and calibration
//-----------------------------------------

var c1, c2, y, R, r, pi, tau, f_tau, b, utility, EV, z;
varexo e, e_a, e_b;

parameters alfa, betta, chi, deltta, dummy_IT, dummy_level, eps, gama, n, pistar, ybar, gbar, bstar, c1star, c2star, utilitystar, taustar, Rstar, phi, sig_e, sig_A;

alfa = 0.3; 
betta = 0.85;  
chi = 1;
deltta = 0.3;
dummy_IT = 1;
dummy_level = 0;  //1 for NGDP-level taregting
gama = 5;   
eps = 0.5;  
n = 0.4;  
pistar = 1.8;
ybar = 1;
gbar = 0.15;  
phi = 0.5;
sig_e = 0.025;  
sig_A = 0.05;  
sig_b = 0.00;  //0.01

//----------------------------------
//1. Find steady state init vals
//----------------------------------
NGDP_steady_state

bstar = b_root;
Rstar = pistar*chi*(1+n);
taustar = ( gbar + (chi-1)*bstar ) / ybar;
c1star = (1-alfa-taustar)*ybar - phi*taustar^2 - bstar;
c2star = alfa*(1+n)*ybar + chi*(1+n)*bstar;
utilitystar = (1/(1-gama))*(c1star^eps + betta*c2star^eps)^((1-gama)/eps);

//--------------------------------
//2. Model
//--------------------------------
model;

//Output
y = ybar*exp(e_a);

//Consumption when young 
c1 = (1-alfa-tau)*y - phi*tau^2 - b;

//Consumption when old
c2 = alfa*(1+n)*y + r*b(-1);

//Bond supply (alternative rules)
//b = bstar;
//b = bstar + deltta*(b(-1) - bstar) + e_b;
//b = bstar - deltta*(y(-1) - ybar) + 0.00000001*e_b;
b = bstar - deltta*(tau(-1)*y(-1) - steady_state(tau)*ybar) + 0.00000001*e_b;

//Determination of taxes 
tau = ( gbar + r*b(-1)/(1+n) - b ) / y;

//Consumption Euler equation 
1 = betta*(R/pi(+1))*(c1/c2(+1))^(1-eps)*( c2(+1)/( EV^(1/(1-gama)) ) )^(1-gama-eps);

//Lagged innovation
z = e(-1);

//Determination of inflation (IT if dummy_IT=1, NGDP if dummy_IT=0)
pi = dummy_IT*pistar*exp(e) + (1-dummy_IT)*pistar*(y(-1)/y)*exp(e-dummy_level*e(-1));

//Real interest rate on bonds 
r = R(-1)/pi;

//Tax burden
f_tau = phi*tau^2;

//Lifetime utility 
utility = 1/(1-gama)*( c1^eps + betta*( EV )^(eps/(1-gama)) )^((1-gama)/eps);

//Expectation term
EV = c2(+1)^(1-gama);

end;

//----------------------------------------
//3. Initial values and shock calibration
//----------------------------------------

initval;
c1 = c1star;
c2 = c2star;
b = bstar;
R = Rstar;
pi = pistar;
y = ybar;
tau = taustar;
utility = utilitystar;
EV = c2^(1-gama);
r = R/pistar;
z = 0;
end;

steady;

shocks;
var e; stderr sig_e;
var e_a; stderr sig_A;
var e_b; stderr sig_b;
end;

//---------------------------
//Find optimal bond supply
//---------------------------

n_loop = 100;
chi_stack = linspace(0.9,1.045,n_loop);
Stack_utility = NaN(n_loop,1); Stack_tau = Stack_utility; Stack_b = Stack_utility; 

for i=1:n_loop

    //----------------------------------------------
    //Find determinstic SS and use in Dynare solver
    //----------------------------------------------
    chi = chi_stack(i);
    NGDP_steady_state

    bstar = b_root;
    Rstar = pistar*chi*(1+n);
    taustar = ( gbar + (chi-1)*bstar ) / ybar;
    c1star = (1-alfa-taustar)*ybar - phi*taustar^2 - bstar;
    c2star = alfa*(1+n)*ybar + chi*(1+n)*bstar;
    utilitystar = (1/(1-gama))*(c1star^eps + betta*c2star^eps)^((1-gama)/eps);

    //---------------------------
    //Stochastic simulations
    //---------------------------
    steady;

    stoch_simul(order=2, drop=0, periods=0, irf=0, noprint);
    
    Loop_record  //Record results in each loop

end

[U_ss_max,Index_ss] = max(U_ss);
[U_stoch_max,Index_stoch] = max(Stack_utility);

NGDP_Baseline_Plotter

//--------------------
//Optimal bond policy
//--------------------
b_ss_opt = b_ss(Index_ss)
r_ss_opt = r_ss(Index_ss)

b_opt = Stack_b(Index_stoch);
r_opt = Stack_r(Index_stoch);

//------------------------
//Welfare at the optimum
//------------------------
if dummy_IT == 1
    Utility_max_IT = Stack_utility(Index_stoch);
else 
    Utility_max_NGDP = Stack_utility(Index_stoch);  
    Lambda_NGDP = 100*( (Utility_max_NGDP/Utility_max_IT)^(1/(1-gama)) - 1)  //CE welfare gain under NGDP
end








