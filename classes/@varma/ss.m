function SSR = ss(obj,V,sid)
% Function SS
%
% Purpose:    Represent VARMA process x(t) in state space form
%
% Format:     SSR = ss(obj,V,sid)
%
% Input:      obj       varma object
%             V         innovation cov matrix
%             sid       signal ID of x(t)
%
% Output:     SSR       state space representation (structure)
%
%                       y(t) = SSR.F*y(t-1) + v(t), E(vv') = SSR.Q
%                       x(t) = SSR.H*y(t)
%
%                       SSR.P - s.s. MSE cov matrix for predicting y
%
%                       P = F[P - P*H'*inv(H*P*H')*H*P]*F' + Q
%
%                       SSR.K - s.s. Kalman gain
%
%                       K = F*P*H'*inv(H*P*H')
%
%                       SSR.V - Wold innovation cov matrix
%
%                       V = H*P*H'
%
% Written by Fei Tan, Saint Louis University
% Updated: September 30, 2020

%% -------------------------------------------
%          VARMA to State Space Form
%---------------------------------------------

% Initialization
p = length(obj.AR);
q = length(obj.MA);
m = size(obj.MA{1},1);
ns = length(sid);

% Construct state space form
r = max(p,q);
AR = cell2mat(obj.AR(:));
AR(end+1:r*m,:) = 0;
MA = cell2mat(obj.MA(:));
MA(end+1:r*m,:) = 0;
SSR.F = [AR [eye(m*(r-1));zeros(m,m*(r-1))]];
SSR.Q = MA*V*MA';
SSR.H = zeros(ns,m*r);
for k = 1:ns
    SSR.H(k,sid(k)) = 1;
end

% Solve discrete-time Riccati
A = SSR.F';
B = SSR.H';
U = cholmod(SSR.Q);
Q = U'*U;
R = zeros(ns);
S = zeros(m*r,ns);
E = eye(m*r);
[SSR.P,K,~] = idare(A,B,Q,R,S,E);
SSR.K = K';
SSR.V = SSR.H*SSR.P*SSR.H';

%-------------------- END --------------------
