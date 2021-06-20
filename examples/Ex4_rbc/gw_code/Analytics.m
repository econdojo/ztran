H_xi_rk=lambda3*(nu*(1-alpha)-1);
H_xi_ra= lambda3*(1+alpha*nu);
H_c_r = - lambda3*nu;

H_xi_wk = (1-alpha)*(1 - (1-alpha)*nu);
H_xi_wa = alpha*(1 - (1-alpha)*nu);
H_c_w = (1-alpha)*nu;

H_xi_kk = 1;
H_xi_theory = [ H_xi_rk     H_xi_ra     0     0
                H_xi_wk     H_xi_wa     0     1
                1           0           1     0]';
 
H_c_theory = [H_c_r H_c_w   0]';

H_base';
H_xi_theory;
H_c_theory;

F_xi_kk = lambda1 + lambda2*nu*(1-alpha);
F_xi_ka = lambda2*(1+alpha*nu);
F_c_k = lambda4 - nu*lambda2;

F_xi_ksks = lambda1;
F_xi_ksas = lambda2*(1+squig);
F_c_ks = lambda4 - squig*lambda2;


F_xi_theory = [ F_xi_kk     0           F_xi_ka     0
                0           F_xi_ksks   0           F_xi_ksas];


F_c_theory = [F_c_k -F_c_ks]';
F_s_theory = [0 F_c_ks]';

%F_xi_theory
%F_c_theory
 


