%Loop record
%Insert in NGDP_2024_Baseline.mod and its variants

    Stack_c1(i) = oo_.mean(1); 
    Stack_c2(i) = oo_.mean(2);
    var_c1(i) = oo_.var(1,1);
    var_c2(i) = oo_.var(2,2);
    Stack_tau(i) = oo_.mean(7);
    var_tau(i) = oo_.var(7,7);
    Stack_ftau(i) = oo_.mean(8);

    Stack_b(i) = oo_.mean(9);
    Stack_r(i) = oo_.mean(5);
    var_r(i) = oo_.var(5,5);
    Stack_utility(i) = oo_.mean(10);

    b_ss(i) = oo_.steady_state(9);
    r_ss(i) = oo_.steady_state(5);
    U_ss(i) = oo_.steady_state(10);
    Lambda_ss(i) = 100*( (U_ss(i)/U_ss(1))^(1/(1-gama))  - 1);   %Consumption equiv. welfare gain 

