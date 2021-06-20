etadiff = 1;
EtaIts = 0;

while etadiff>ITERATE_ACCURACY
    
    last_eta = eta;
    
    eta = CalcEta(eta, mu_w,mu_c, L,M);

    if RESTRICTED_ETA
        eta = RestrictEta(eta);
    end

    etadiff = max(max(abs(eta-last_eta)));    
    
    EtaIts = EtaIts + 1;
    if EtaIts > MAX_ITS
        error ('eta failed to converge')
    end
    
    if It_Tick(EtaIts)
        text = ['Iteration: ' num2str(AllIts) ', Eta Iteration '  num2str(EtaIts) '  ' num2str(etadiff) ];
        disp(text);
    end
end

if RESTRICTED_ETA
    eta = RestrictEta(eta);
end

