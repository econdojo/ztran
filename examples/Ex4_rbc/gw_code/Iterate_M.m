Mdiff = 1;
MIts = 0;

%M_Term = I_states-beta * H;

while Mdiff >ITERATE_ACCURACY
    last_M = M;
    
    M = CalcM(L,eta,beta,H,M,MW);

    Mdiff  = max22(M,last_M);  
    MIts = MIts + 1;
    
    if MIts > MAX_ITS
        M
        error ('M failed to converge')
    end
end

if DEBUG
    MIts
end
