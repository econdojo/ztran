% iterate beta and M together

MBetaEtadiff = 1;
MBetaEtaIts = 0;
while MBetaEtadiff >ITERATE_ACCURACY
    last_M = M;
    last_beta = beta;
    last_eta = eta;
    
    H = CalcH(H_base, H_c, H_s,eta);
    
    [beta P] = CalcBeta(P,H,M,N*Q*N',beta);

    MW = [F_w F_c*eta'];
    M = CalcM(L,eta,beta,H,M,MW);
    
    N = CalcN(beta,H,N);

    eta = CalcEta(eta, mu_w,mu_c, L,M);

    if RESTRICTED_ETA
        eta = RestrictEta(eta);
    end
    
    MBetaEtadiff =max22(beta,last_beta)+max22(M,last_M)+max22(eta,last_eta); 
    
    MBetaEtaIts = MBetaEtaIts + 1;
    
    if It_Tick(MBetaEtaIts)
        text = ['M,Beta, eta Iteration '  num2str(MBetaEtaIts) '  ' num2str(MBetaEtadiff) ];
        disp(text);
    end

   
    if MBetaEtaIts > MAX_ITS
        error ('M and Beta iteration failed to converge')
    end
end

if DEBUG
    MBetaIts
end
