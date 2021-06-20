global eta_full_inf

diff = 1;
Its = 0;
eta_f = ones(r,2);
MaxEtaIts=100000;
F_c_eta = [F_c+F_s F_s];

eta_f = ones(r,2);
    
while diff>ITERATE_ACCURACY/10000
    last_eta_f = eta_f;
    
    mu_eta = [mu_w+ [mu_c mu_s]  * eta_f'
             zeros(1,r)];

    eta_f = ((eta_f'-sigma_c*mu_eta)*(F_w+ F_c_eta*eta_f'))';
    maxdiff = 0;
    diff = max(max(abs(eta_f-last_eta_f)));     
    Its = Its + 1;
    if Its > MaxEtaIts
        error ('eta failed to converge')
    end
end


eta_full_inf = eta_f(:,1) + eta_f(:,2);

