% iterate beta and M together

MBetadiff = 1;
MBetaIts = 0;


while MBetadiff >ITERATE_ACCURACY
    last_M = M;
    last_beta = beta;

    M = CalcM(L,eta,beta,H,M,MW);
    N = CalcN(beta,H,N);

    [beta P] = CalcBeta(P,H,M,N*Q*N',beta);

    MBetadiff =max22(beta,last_beta)+max22(M,last_M);
    
    MBetaIts = MBetaIts + 1;
    
    if It_Tick(MBetaIts)
        text = ['Iteration: ' num2str(AllIts) ', M and Beta Iteration '  num2str(MBetaIts) '  ' num2str(MBetadiff) ];
        disp(text);
    end
   
    if MBetaIts > MAX_ITS
        error ('M and Beta iteration failed to converge')
    end
end

if DEBUG
    MBetaIts
end
