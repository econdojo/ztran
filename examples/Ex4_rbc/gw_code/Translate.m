% take system in form
% W(t+1) = G_w W(t) + G_c c(t) + G_s cs(t) + G_y ys(t)
% M_y y(t) = M_w W(t) + M_c c(t) + M_s cs(t)

% and translate it into form
% W(t+1) = F_w W(t) + F_c c(t) + F_s cs(t) + G_y ys(t)
% i(t) = H_w W(t) + H_c c(t) + H_s cs(t)

global I_states n r

n_y = size(M_y,2);
n= size(OBSERVABLES,2); 
r= size(G_w,1); 
r_xi = size(Xi_agg,2);
r_chi = r - r_xi;

F_w = G_w + G_y*inv(M_y) * M_w;
F_c = G_c + G_y*inv(M_y) * M_c;
F_s = G_s + G_y*inv(M_y) * M_s;


% invert
H_w_all =   inv(M_y) * M_w;
H_c_all =   inv(M_y) * M_c;
H_s_all =  inv(M_y) * M_s;
H_om_all =  inv(M_y) * M_om;

% pick out rows corresponding to OBSERVABLES
n = size(OBSERVABLES,2);

H_w = H_w_all(OBSERVABLES,:);
H_c = H_c_all(OBSERVABLES,:);
H_s = H_s_all(OBSERVABLES,:);
H_om = H_om_all(OBSERVABLES,OBSERVABLES);


mu_w = zeros (1,r);

mu_w(1,:)= H_w_all(r_index,:);
mu_c = H_c_all(r_index);
mu_s = H_s_all(r_index);

Q = Q_w;

% add processes for measurment error under sthe statez
for i = 1:size(MEAS_ERROR,2)
        ObsTemp = MEAS_ERROR(i);
        F_w = [   F_w             zeros(r,1) 
                  zeros(1,r)   zeros(1,1)]   ;
              
        F_c = [F_c' 0]';
        F_s = [F_s' 0]';
        
        for j = 1:size(OBSERVABLES,2)
            if OBSERVABLES(j) == ObsTemp
                break
            end
        end
        H_w = [ H_w  H_om(:,j)];
        H_w_all = [ H_w_all zeros(n_y,1)];

        mu_w = [mu_w 0];
        
        Q = [   Q            zeros(r,1)
                zeros(1,r)   var_meas];
    
        r = r + 1;               % add measurement erros onto states
end

L = [   F_s
        zeros(h*r,1)];
    
H_base = [H_w zeros(n, h*r)];

I_states = eye((h+1)*r);

