function obj = fit(z,fz,p,q,or)
% Function FIT
%
% Purpose:    Fit joint VARMA representation to data points (note: fit separate 
%             VARMA's for aggregate & individual variables to improve accuracy)
%
% Format:     obj = varma.fit(z,fz,p,q)
%
% Input:      z(k)      point of evaluation
%             fz(:,:,k) VARMA representation evaluated at z(k)
%             p         VAR order
%             q         VMA order
%             or        order reduction for stationarity (logical)
%
% Output:     obj       varma object
%
% Written by Fei Tan, Saint Louis University
% Updated: September 30, 2020

%% -------------------------------------------
%                  Fit VARMA
%---------------------------------------------

% G0*x = G1
z = z(:);
[m,n,N] = size(fz);
G1 = reshape(fz,m,n*N)';

if p==0
    AR = [];
else
    AR = zeros(n*N,m*p);
    zz = kron(z,ones(n,m));
    AR(:,1:m) = zz.*G1;
    for k = 2:p
        AR(:,((k-1)*m+1):(k*m)) = zz.*AR(:,((k-2)*m+1):((k-1)*m));
    end
end

MA = zeros(n*N,n*(q+1));
zz = kron(z,ones(n));
MA(:,1:n) = repmat(eye(n),N,1);
for k = 2:(q+1)
    MA(:,((k-1)*n+1):(k*n)) = zz.*MA(:,((k-2)*n+1):((k-1)*n));
end

% VARMA representation
while true
    G0 = [AR MA];
    [U,S,V] = svdc(G0);
    x = (V/S*U'*G1)';
    
    if p>1 && or
        A = [x(:,1:(p*m));eye((p-1)*m) zeros((p-1)*m,m)];
    else
        break
    end
    
    if any(abs(eig(A))>1)    % stationarity test
        p = p-1;
        if q>0
            q = q-1;
        end
        AR = AR(:,1:(p*m));
        MA = MA(:,1:((q+1)*n));
    else
        break
    end
end

AR = reshape(x(:,1:(p*m)),m,m,p);
if isempty(AR)
    AR = {};
else
    AR = squeeze(num2cell(AR,[1 2]));
end

MA = reshape(x(:,(p*m+1):end),m,n,q+1);
if isempty(MA)
    MA = {};
else
    MA = squeeze(num2cell(MA,[1 2]));
end
obj = varma(AR,MA);

%-------------------- END --------------------
