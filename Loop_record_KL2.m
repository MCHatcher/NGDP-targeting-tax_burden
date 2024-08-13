%Loop_record_KL2
%Used in the file NGDP_2024_KL2_lump_sum.mod and its variants

    %Risky steady state

    Stack_c1(j) = oo_.mean(1); 
    Stack_c2(j) = oo_.mean(2);
    var_c1(j) = oo_.var(1,1);
    var_c2(j) = oo_.var(2,2);
    Stack_tau(j) = oo_.mean(12);
    var_tau(j) = oo_.var(12,12);
    Stack_comp(j) = oo_.mean(29);

    Stack_y(j) = oo_.mean(5);
    Stack_k(j) = oo_.mean(3);
    Stack_l(j) = oo_.mean(6);
    Stack_rk(j) = oo_.mean(9);
    var_k(j) = oo_.var(3,3);
    var_rk(j) = oo_.var(9,9);

    Stack_b(j) = oo_.mean(4);
    Stack_r(j) = oo_.mean(8);
    var_r(j) = oo_.var(8,8);
    Stack_utility(j) = oo_.mean(13);

    %Stack_R(j) = oo_.mean(7);

    %With LS taxes
    Stack_utility1(j) = oo_.mean(25);

    %Deterministic steady state

    Stack_bss(j) = oo_.steady_state(4);
    Stack_rss(j) = oo_.steady_state(8);
    Stack_Uss(j) = oo_.steady_state(13);
    Stack_rss1(j) = oo_.steady_state(21);
    Stack_Uss1(j) = oo_.steady_state(25);
    
    Stack_Lambda_ss(j) = 100*( (Stack_Uss(j)/Stack_Uss1(1))^(1/ (thetta*(1-gama)) )  - 1);   %Consumption equiv. welfare gain 
    Stack_Lambda_ss1(j) = 100*( (Stack_Uss1(j)/Stack_Uss1(1))^(1/ (thetta*(1-gama)) )  - 1);   %Consumption equiv. welfare gain
    Tax_burden_ss(j) = 100*( (Stack_Uss1(j)/Stack_Uss(j) )^(1/( thetta*(1-gama) ) ) -1);

    if dummy_IT==1
        U_IT_mean(j) = oo_.mean(13);
        Tax_burden_welfare_loss(j) = 100*( (Stack_utility1(j)/Stack_utility(j) )^(1/( thetta*(1-gama) ) ) -1);
    else
        Lambda_NIT(j) = 100*( (U_IT_mean(j)/Stack_utility(j))^(1/( thetta*(1-gama) ) ) -1);
        Tax_burden_welfare_loss_NIT(j) = 100*( (Stack_utility1(j)/Stack_utility(j) )^(1/( thetta*(1-gama) ) ) -1);
    end


