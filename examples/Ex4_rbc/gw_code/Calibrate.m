% calcuates the properties of the process for idiosyncratic labour income
% regression to find persistence of process

% run NumRegress regressions of PERIODS length each then caluclatyes
% average annual AR coefficient and innocvation standar ddevaition

NumRegress = 100;

SIMULATE = 1;
PERIODS = 10000;

Sample =[100 200 500 1000 10000 ];
if max(Sample)>PERIODS
    error('Max sample');
end


RHO = zeros(size(Sample,2),NumRegress);
INNOV = zeros(size(Sample,2),NumRegress);

for Regress = 1:NumRegress
    % generate time series
    Simul
    
    % labour is in deviations so need to add on aggregate
    income = Controls(i_lsd,:) + Controls(i_wsd,:);
    
    for i = 1:size(Sample,2)
        PERIODS = Sample(i);
        X = income(3:PERIODS);
        X_lag = income(2:PERIODS-1);
    
        % regress income on its lag
        rho_q = inv(X_lag*X_lag')*X_lag*X';

        % calcute innovation standard deviation and compare it with Guvenen's estimate
        innov_q = zeros(1,PERIODS);
        for regress_i=2:PERIODS
            innov_q(regress_i) = income(regress_i) - rho_q*income(regress_i-1);
        end

        sd_q = std(innov_q);    
    
        % calcuate and save annual values
        INNOV(i,Regress) = sd_q*(1+rho_q^2+rho_q^4+rho_q^6)^0.5;
        RHO(i,Regress) = rho_q^4;
    
    end
    if (Regress/10 == int32(Regress/10));
        Regress
    end
end

format short
disp('Average annual ar parameter')
disp(Sample)
disp(mean(RHO'))
disp('Average annual innovation standard deviation')
disp(Sample)
disp(mean(INNOV'))

% check against Guvenen's estimate
disp('Guvenen estimate of annual standard deviation')
disp(0.13/0.007)
