function [Z_sigma] = CalcUncondVarIter(A, B,Q)

% calculates the unconditional variance of a system
% Z(t) = A Z(t-1) + B v(t)
% where A is a matrix, v is a vector of shocks 

max_its = 10000;
its = 0;
accuracy = 10e-5;
Z_sigma =zeros(size(A));
diff = 1;

while diff > accuracy
    last = Z_sigma;
    Z_sigma = A*last*A' + B*Q*B';
    % zero out explosive element
    Z_sigma(3,3) = 0;
    diff = max(max(abs(Z_sigma-last))) - abs(Z_sigma(3,3) - last(3,3));
    its = its +1;
    if its>max_its
        error('Unconditional vairance iteration failed to converge');
    end
end

% recalc with explosive element
Z_sigma = A*last*A' + B*Q*B';







