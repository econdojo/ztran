% produces simulated time series 
% expects matrices for VAR:  Z(t) = AR Z(t-1) + NN nu(t)


%****************************************************************************************************************  
%****************************SET UP VARIABLES                   *************************************************  
%****************************************************************************************************************  
C_agg = zeros(1,PERIODS);
C_certain = zeros(1,PERIODS);
C_idio = zeros(1,PERIODS);
N_agg = zeros(1,PERIODS);
N_idio = zeros(1,PERIODS);
W_idio = zeros(1,PERIODS);
Z_s = zeros ((h+1)*r+h*r,PERIODS);
Z_full = zeros (r,PERIODS);
C_full = zeros(1,PERIODS);;
Nu = zeros (r,PERIODS);
Controls = zeros(n_y,PERIODS);

if SIMULATE == 1
    Nu(2,:) = randn(1,PERIODS)*sqrt(var_agg);
    Nu(4,:) = randn(1,PERIODS)*sqrt(var_idio);
else
    % aggregate tech shock
    Nu(2,1) = 1;
    % idio tec shock
    Nu(4,1) = 0;
end

%****************************************************************************************************************  
%**************************** TIME SERIES                *************************************************  
%****************************************************************************************************************  


for i = 2:PERIODS
    Z_s(:,i) = AR*Z_s(:,i-1) + NN*Nu(:,i-1);       
    C_agg(i) = eta'*T2*Z_s(1:(h+1)*r,i);
    C_idio(i) = eta'*Z_s((h+1)*r+1:(h+1)*r+h*r,i);
    Controls(:,i) = H_w_all*Z_s(1:r,i)+H_c_all*C_agg(i)+H_s_all*C_idio(i);
    if PLOT_CERTAIN | PLOT_CERTAIN_DEV | HIERARCHY_IMPACT
        C_certain(i) = eta_full_inf(1:r)'*Z_s(r+1:2*r,i);
    end
    if PLOT_FULL_INF
           Z_full(:,i) = (F_w + (F_s +F_c) * [eta_full_inf(1:2)' zeros(1,r-2)]+ F_s * [0 0 eta_full_inf(3:r)'])* Z_full(:,i-1)+ Nu(1:r,i-1);
           C_full(i) = [eta_full_inf(1:2)' zeros(1,r-2)]*Z_full(:,i);
           Controls_full(:,i) = H_w_all*Z_full(:,i)+H_c_all*C_full(i);
    end   
    
    %if (i/1000 == int32(i/1000))
    %    i
    %end
end


if Z_s(3,PERIODS) > ITERATE_ACCURACY & Nu(4,1) == 0 
    disp(['Long run idiosync capital ' num2str(Z_s(3,PERIODS))])
    pause   
end



            




