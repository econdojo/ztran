
betadiff = 1;
BetaIts = 0;

NQN = N*Q*N';

while betadiff>ITERATE_ACCURACY
    last_beta = beta;
    last_P = P;
    
    [beta P] = CalcBeta(P,H,M,NQN,beta);
     
    betadiff =max22(beta,last_beta);
    BetaIts = BetaIts + 1;
    
    if It_Tick(BetaIts)
        text = ['Beta Iteration ' num2str(AllIts) '  ' num2str(BetaIts) '  ' num2str(betadiff) ];
        disp(text);
    end
    
    if BetaIts > MAX_ITS
        error('beta failed to converge');
    end
end

