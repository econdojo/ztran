% takes model in form

% defined as global so can be used in functions 
global r h T1 T2 T3 T7 T8 T11 T12
% defined as global so can be used in functions from the command line
global eta Z_s

%***********************************************************************************************************
%****************************** truncation matrices      ******************************************************
%***********************************************************************************************************%
T1 = [  eye(h*r)
        zeros(r,(h-1)*r) eye(r)];
    
T2 = [  zeros(h*r,r)       eye(h*r)];

T3 = [eye(h*r) zeros(h*r,r)];


% sets idiosyncratic elements of W to zero
% if measurement error is in an idiosyncratic variable it is averaged
% out by T8
T8 = zeros(r);
T8(1:r_xi,1:r_xi) = eye(2);
for i = 1:size(MEAS_ERROR)
    temp = MEAS_ERROR(i);
    for j=Y_agg
        if temp == j
            T8(r_xi+r_chi+i,r_xi+r_chi+i) =1;
        end
    end
end

% sets idiosyncratic elements of X to zero - not hierachy is all aggregate
    
T7 = [  T8              zeros(r,h*r)
        zeros(h*r,r)    eye(h*r) ];
    
%***********************************************************************************************************
%******************************starting values        ******************************************************
%***********************************************************************************************************
Iterate_FullInfEta

if LOAD_VARS
    load('StartVars'); 
    % keep copies of the starting values
    M_start = M;
    P_start = P ;
    beta_start = beta ;
    eta_start = eta  ;
    N_start =N  ;
else
    StartValues 
end

M = M_start;
P = P_start;
beta = beta_start;
eta = eta_start;
N = N_start;    
%***********************************************************************************************************
%****************************** restrict state vector    ******************************************************
%***********************************************************************************************************

if RESTRICTED_ETA == 1 %& max(max(H_c))== 0
    % fixed labour case
    %if max(max(H_c))~= 0
    %    error('Eta restrictions do not work for endogenous labour');
    %end
    
    if MEAS_ERROR ~= [i_ws]
        error('Eta restrictions only implemented for mesurement error in wage');
    end

    H_w_temp = H_w(:,1:r_xi+r_chi);
    H1 = [  H_w_temp
            zeros(1,r_xi+r_chi)];

    H1(n+1,2) = 1;

    H2 = [  H_w_temp  zeros(n,1)
            zeros(1,r_xi+r_chi+1)];
        
    H2(n+1,r_xi+r_chi+1) = 1;
    
    % maps [xi a] onto W
    T_xi = zeros(5,3);
    T_xi(1,1) = 1;
    T_xi(2,2) = 1;
    T_xi(5,3) = 1;
    
    % maps W_s onto Xi
    T_w = zeros(2,5);
    T_w(1,1) = 1;
    T_w(2,2) = 1;
    
    % addson row to W for measurtement error
    T_om = zeros(5,4);
    T_om(1:4,1:4) = eye(4);
    
    H_31 = T_om*inv(H1)*H2*T_xi ;
    
    H_41 = H_31(:,1:2);
    H_51 = H_31(:,3);
    
    T11 = [eye(r)   zeros(r,h)];
    NewRows =  [ H_41*T_w H_51  zeros(r,h-1) ];
    T11 = [ T11
            NewRows];
    
    H_3_last = H_31;

    for i = 1:h-1
        
        Temp = [    T_w*H_3_last                        zeros(size(T_w*H_3_last,1),1)];
        for j = 1:i
              NewRow = [zeros(1,size(T_w*H_3_last,2)+1)];
              NewRow(size(T_w*H_3_last,2)-i+1+j) = 1;
              Temp = [  Temp
                        NewRow];
        end
        
        H_3_next  = H_3_last*Temp;
        H_4 = H_3_next(:,1:2);
        H_5 = H_3_next(:,3:i+3);
        NewRows = [ H_4*T_w H_5 zeros(r,h-1-i) ];
        T11 = [ T11
                NewRows];

        H_3_last = H_3_next;
    end
    
    % T12 maps contracts full state vector to reduced state vector
    T_temp = zeros(5,1);
    T_temp(2) = 1;
    T12 = [   eye(r)              zeros(r,h)];
    NewRow = [zeros(r,r)        T_temp         zeros(r,(h-1))];
    T12 = [ T12
            NewRow];
        
    for i = 1:h-1
            NewRow = [zeros(r,r)    zeros(r,i)  T_temp  zeros(r,(h-1-i))];
            T12 = [ T12
                    NewRow];
    end
end
            


%***********************************************************************************************************
%**************  Iteration                                        ******************************************
%***********************************************************************************************************

% IT_STYLE = 1: iterate M until consistent with itself, iterate beta until
% consistent with itself and M, iterate eta until consistent with given M and
% beta, then loop
% IT_STYLE = 2: calculate M and beta, and iterate these two until
% consistent, then calculate eta and loop
% IT_STYLE = 3: calcualte M, beta and eta.  Iterate until mutually
% consistent
%IT_STYLe = 4: iterate M,beta,eta.  Iterate until all consistent.

disp('Iteration starting...');

if IT_STYLE == 3
    Iterate_MBetaEta
elseif IT_STYLE == 4
    Diff4 = 1;
    AllIts = 0;

    MW = [F_w F_c*eta'];
    H = CalcH(H_base, H_c, H_s,eta);

    while Diff4 >ITERATE_ACCURACY
        last_4_M = M;
        last_4_beta = beta;
        last_4_eta = eta;
        MW = [F_w F_c*eta'];
        H = CalcH(H_base, H_c, H_s,eta);
        Iterate_Beta
        Iterate_M
        N = CalcN(beta,H,N);
        Iterate_Eta

        Diff4  = max22(M,last_4_M)+max22(beta,last_4_beta)+max22(eta,last_4_eta); 
        AllIts = AllIts + 1;
        if It_Tick(AllIts) 
            text = ['New Iteration '  num2str(AllIts) '  ' num2str(Diff4) ];
            disp(text);
        end

        if AllIts > MAX_ITS
            error ('New failed to converge')
        end
   end
else
    AllIts = 0;
    AllDiff = 1;
    while AllDiff >ITERATE_ACCURACY
        last_all_eta = eta;
        last_all_beta = beta;
        last_all_M = M;
        last_all_N = N;
    
        MBetaIts = 0;
        MBetaDiff = 1;
    
        MW = [F_w F_c*eta'];
        H = CalcH(H_base, H_c, H_s,eta);

        if IT_STYLE  == 2
            Iterate_MBeta
        elseif IT_STYLE == 1
            % iterate beta and N separately, then iterate until the two
            % converge
            while MBetaDiff >ITERATE_ACCURACY
                last_betaM_M = M;
                last_betaM_N = N;
                last_betaM_beta = beta;

                Iterate_M
                N = CalcN(beta,H,N);
                Iterate_Beta

                MBetaDiff  = max22(M,last_betaM_M)+max22(beta,last_betaM_beta)+max22(N,last_betaM_N); 
                MBetaIts = MBetaIts + 1;
                if It_Tick(MBetaIts) 
                    text = ['M,Beta, eta Iteration '  num2str(MBetaIts) '  ' num2str(MBetaDiff) ];
                    disp(text);
                end

                if MBetaIts > MAX_ITS
                    error ('M / beta failed to converge')
                end
            end
        end
        
        if h>1
            Iterate_Eta
        end
    
        AllDiff  = max22(M,last_all_M)+max22(N,last_all_N)+max22(beta,last_all_beta)+max22(eta,last_all_eta);
        AllIts = AllIts + 1;
        
        if It_Tick(AllIts)
            Text = ['Iteration: ' num2str(AllIts) ', difference: ' num2str(AllDiff)];
            disp(Text);
        end
        if AllIts > MAX_ITS
            error ('Outside iteration to converge')
        end
    end
end


%***********************************************************************************************************
%***************************    Housekeeping                      ******************************************
%***********************************************************************************************************


% check for adding up constraints
sum_eta=SumEta(eta);

eta_diff = sum_eta(1:4) - eta_full_inf(1:4);
disp(' ')
disp(' ')
disp('Difference of sum from full information eta');
disp(eta_diff)


if SAVE_VARS
    save('StartVars','beta','P','M','N','eta');
end


    

%***********************************************************************************************************
%**************  VAR representation                                 ******************************************
%***********************************************************************************************************

% expects matrices for VAR:  Z(t) = AR Z(t-1) + NN nu(t)

AR1 = [  F_w                F_c*eta'                                                    F_s*eta'
        zeros(h*r,r)        T3*(L*eta'+(eye((h+1)*r)-beta * H)*M*T1)                    zeros(h*r,h*r)
        zeros(h*r,r)        zeros(h*r,h*r)                                              T3*(L*eta'+(eye((h+1)*r)-beta * H)*M*T1)];

AR2 = [ zeros(r,(h+1)*r)    zeros(r,h*r) 
        T3*beta*H*M*T7      zeros(h*r,h*r)
        T3*beta*H*M         zeros(h*r,h*r)];
    
AR = AR1+AR2;

NN = [   N
         N(r+1:(h+1)*r,:)];











        
    
