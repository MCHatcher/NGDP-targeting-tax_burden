% Discretization, written by M. Hatcher (m.c.hatcher@soton.ac.uk)
% This code shows how to approximate an IID-Normal process using the
% method in Adda and Cooper (2003, Dynamic Economics) as in the notes 
% of Makoto Nakajima (https://makotonakajima.github.io/comp/).

mu = 0; %Mean of process

%Equal probabilities by construction

%Housekeeping
m_i = NaN(n_states,1); e_i = m_i; 

%Find abscissas and states
for i=1:n_states-1
    m_i(i) = norminv(i/n_states)*sigma + mu;
    
    if i==1
        e_i(i) = mu - sigma*n_states*normpdf( (m_i(i)-mu)/sigma );
    else
        e_i(i) = mu - sigma*n_states*( normpdf( (m_i(i)-mu)/sigma ) -  normpdf( (m_i(i-1)-mu)/sigma ) );
    end
end

e_i(n_states) = mu + sigma*n_states*normpdf( (m_i(n_states-1)-mu)/sigma );


        
        

    




