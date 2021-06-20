function eta = RestrictEta(last_eta)

global h r T3 T11 T12 eta_full_inf FUDGE_MU

mu = (last_eta'*T3*T11)';

if FUDGE_MU
    for i = h-1:h+1
        mu(i) = 0;
    end
end

eta = (T12*mu);

% impose adding up constraint
eta(2) = eta_full_inf(2) - sum(mu(r+1:r+h));


% translate
eta = eta(1:h*r);


