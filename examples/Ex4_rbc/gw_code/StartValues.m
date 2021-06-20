eta_start = [   eta_full_inf
                zeros((h-1)*r,1)];


%***********************************************************************************************************
%******************************starting values for beta ******************************************************
%***********************************************************************************************************
P = eye(r);
diff = 1;
Its = 0;
MaxBetaIts = 100000;

beta = zeros(r,n);
while diff>ITERATE_ACCURACY
    last_beta = beta;
    last_P = P;
    beta = P * H_w'*inv((H_w*P*H_w'));
    MM = (eye(r)-beta * H_w)*P;
    P = F_w * MM * F_w' + Q;
        
    diff_beta = max(max(abs(beta-last_beta)));
    diff_P = max(max(abs(P-last_P)));
    diff = max(diff_P,diff_beta);
    
    Its = Its + 1;
    if Its > MaxBetaIts
        errtext = 'beta / J failed to converge';
        if DEBUG
            warning(errtext);
            return
        else
            error (errtext)
        end
    end
end

beta_start = [  beta
                zeros(h*r,n)];;

    
%***********************************************************************************************************
%******************************starting values for M ******************************************************
%***********************************************************************************************************

MW = [F_w F_c*eta_start'];

M_start = [ MW
            zeros(h*r,(h+1)*r)];
        
N_start = [     eye(r)
                T3*beta_start*H_w*T8];
    
P_start = eye((h+1)*r,(h+1)*r);
P_start = zeros((h+1)*r,(h+1)*r);
P_start(1:r,1:r)=P;

M = M_start;
P = P_start;
beta = beta_start;
eta = eta_start;
N = N_start;



      

    
