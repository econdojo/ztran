function eta = CalcEta(last_eta, mu_w,mu_c, L,M)

global r h T1 T2 T3 sigma_c

R = ([mu_w zeros(1,h*r)] + mu_c*last_eta'*T2)';
eta = (((last_eta'*T3 - sigma_c*R')*(L*last_eta'*T3+M))*T1)';



    