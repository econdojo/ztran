function res = irf(obj,imp,crit)
% Function IRF
%
% Purpose:    Simulate VARMA impulse response function
%
% Format:     res = irf(obj,imp)
%
% Input:      obj       varma object
%             imp(:,t)  impulse at period t
%             crit      criterion for nonzero responses
%
% Output:     res(:,t)  response at period t
%
% Written by Fei Tan, Saint Louis University
% Updated: September 30, 2020

%% -------------------------------------------
%           VARMA Impulse Response
%---------------------------------------------

% Initialization
if nargin<3
    crit = 1e-5;
end
p = length(obj.AR);
q = length(obj.MA)-1;
[m,n] = size(obj.MA{1});
N = size(imp,2);
imp = [zeros(n,q) imp];
res = zeros(m,N+p);

% Impulse response
for k = 1:N
    res(:,p+k) = obj.MA{1}*imp(:,q+k);
    for j = 1:p
        res(:,p+k) = res(:,p+k)+obj.AR{j}*res(:,p+k-j);
    end
    for j = 1:q
        res(:,p+k) = res(:,p+k)+obj.MA{j+1}*imp(:,q+k-j);
    end
    
end
res = res(:,p+1:end);
res(abs(res)<crit) = 0;

%-------------------- END --------------------
