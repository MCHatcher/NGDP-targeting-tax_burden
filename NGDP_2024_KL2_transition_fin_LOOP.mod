//Dynare code for simulating transition to NIT in the extended model (Loop version). 
//Many simulations with generational utilities averaged over all simulations. 
//Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). Any errors are my own.

//-----------------------------------------
//1. Variable declaration and calibration
//-----------------------------------------

var c1, c2, k, b, y, l, R, r, rk, w, pi, tau, utility, EV, x_TR, x_c1;
varexo e, e_a;

parameters alfa, betta, chi, dummy_IT, eps, gama, tau_p, thetta, n, pistar, psi, gbar, bstar, c1star, c2star, utilitystar, taustar, Rstar, 
sig_e, sig_A; 

alfa = 0.3; 
betta = 0.85; 
dummy_IT = 1;
gama = 5; 
eps = 0.5;  
n = 0.4;  
pistar = 1.8;
gbar = 0.05;  
thetta = 0.45;  
psi = 1;   
tau_p = 0.025;  
sig_e = 0.025;  
sig_A = 0.05;  

chi = 1.27785714285714;  //Gives (approx.) optimal bond supply; 

//----------------------------------
//1. Find steady state init vals
//----------------------------------
Steady_state_KL2_insert

bstar = b_ss; kstar = k_ss; 
Rstar = pistar*chi;
s_star = s_ss;
taustar = tau_ss;
c1star = c1_ss; c2star = c2_ss;
lstar = l_ss;
utilitystar = utility_ss;

//--------------------------------
//2. Model
//--------------------------------

model;

//Output
y = exp(e_a)*(k(-1)/(1+n))^alfa*l^(1-alfa);

//Consumption when young 
c1 = (1-tau-tau_p)*w*l - k - b;

//Consumption when old
c2 = (1-psi*tau)*rk*k(-1) + R(-1)/pi*b(-1) + x_TR;

//Pension transfer
x_TR = tau_p*(1+n)*w*l;

//Bond supply
b = bstar;

//Determination of taxes 
tau = ( gbar + r*b(-1)/(1+n) - b ) / ( w*l + psi*rk*(k(-1)/(1+n)) );

//Consumption Euler equation (bonds)
1 = betta*(R/pi(+1))*(c2(+1)/c1)^(eps*thetta-1)*( 1/(1-l) )^( (1-thetta)*eps )*( c2(+1)^thetta / EV^(1/(1-gama)) )^(1-gama-eps);

//Consumption Euler equation (capital) 
1 = betta*(1-psi*tau(+1))*rk(+1)*(c2(+1)/c1)^(eps*thetta-1)*( 1/(1-l) )^( (1-thetta)*eps )*( c2(+1)^thetta / EV^(1/(1-gama)) )^(1-gama-eps);

//Labour supply
thetta*(1-l) = (1-thetta)*c1/( (1-tau-tau_p)*w );

//Determination of inflation (IT if dummy_IT=1, NGDP if dummy_IT=0)
pi = dummy_IT*pistar*exp(e) + (1-dummy_IT)*pistar*(y(-1)/y)*exp(e);

//Real interest rate on bonds 
r = R(-1)/pi;

//Return on capital 
rk = alfa*y/(k(-1)/(1+n));

//Wage
w = (1-alfa)*y/l;

//Lifetime utility 
utility = (1/(1-gama))*( ( c1^thetta*(1-l)^(1-thetta) )^eps + betta*( EV )^(eps/(1-gama)) )^((1-gama)/eps);

//Expectation term
EV =  c2(+1)^(thetta*(1-gama));

//Composite consumption 
x_c1 = c1^thetta*(1-l)^(1-thetta);

end;

//----------------------------------------
//3. Initial values and shock calibration
//----------------------------------------

initval;
c1 = c1star;
c2 = c2star;
k = kstar;
b = bstar;
pi = pistar;
R = pistar*chi;
l = lstar;
y = (kstar/(1+n))^alfa*lstar^(1-alfa);
w = (1-alfa)*y/l;
tau = taustar;
utility = utilitystar;
EV = c2star^(thetta*(1-gama));
r = R/pistar;
rk = r/(1-psi*tau);
x_TR = tau_p*(1+n)*w*l;
x_c1 = c1^thetta*(1-l)^(1-thetta);
end;

steady;

histval;
k(0) = 0.0895543096920283;
y(0) = 0.321122804333936;
R(0) = 2.28345163325862;
end;

shocks;
var e; stderr sig_e;
var e_a; stderr sig_A;
end;

xc1_init = 0.187713523038578;
U_init_sum = 0; Utility_sum = zeros(50,1);
n_sim = 50000; 

for j=1:n_sim

    set_dynare_seed(j)
    stoch_simul(order=2, drop=0, periods=50, irf=0, noprint);  
    T = length(oo_.endo_simul(1,1:end));

    c2_init = oo_.endo_simul(2,1);
    U_init = (1/(1-gama))*( ( xc1_init )^eps + betta*( c2_init )^eps )^((1-gama)/eps);
    U_init_sum = U_init_sum + U_init;

    Util = oo_.endo_simul(13,1:end)';
    Utility_sum = Utility_sum + Util;

end

if dummy_IT==1

    U_init_IT = mean(U_init);
    U_mean = Utility_sum/n_sim;
    U_IT = U_mean;

else

    U_init_NIT = mean(U_init);
    U_mean = Utility_sum/n_sim;
    U_NIT = U_mean;

    for t=1:T
        Lambda_NIT(t,1) = 100*( (U_NIT(t)/U_IT(t))^(1/(thetta*(1-gama)) ) - 1);   
    end
    Lambda_init = 100*( (U_init_NIT/U_init_IT)^(1/(thetta*(1-gama)) ) - 1);  
Generations = 0:10;
figure(1) 
hold on, plot(Generations,[Lambda_init; Lambda_NIT(1:10,1)],'k','LineWidth', 1), hold on,
xlabel('Generations (d.o.b.)'), ylabel('% c.e. welfare gain'); yline(0)
end



 