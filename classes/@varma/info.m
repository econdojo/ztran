function [Hu,Hc,P0,P1,P2] = info(obj,V,sid)
% Function INFO
%
% Purpose:    Compute Shannon information flow
%
% Format:     [Hu,Hc,P0,P1,P2] = info(obj,V,sid)
%
% Input:      obj       varma object
%             V         innovation cov matrix
%             sid       signal ID
%
% Output:     Hu        unconditional entropy
%             Hc        conditional entropy
%             P0        unconditional cov matrix
%             P1        s.s. updated MSE cov matrix
%             P2        s.s. predictive MSE cov matrix
%
% Note:       Hu-Hc = Shannon information flow
%
% Written by Fei Tan, Saint Louis University
% Updated: September 30, 2020

%% -------------------------------------------
%                   Entropy
%---------------------------------------------

% Unconditional entropy
P0 = real(integral(@(w)spectrum(w,obj,V),-pi,pi,'ArrayValued',true));% inversion formula
Hu = log2(2*pi*exp(1)*diag(P0))/2;

% Conditional entropy
m = size(obj.MA{1},1);
SSR = ss(obj,V,sid);
P2 = SSR.P(1:m,1:m);
P = SSR.P-(SSR.P*SSR.H')/(SSR.H*SSR.P*SSR.H')*SSR.H*SSR.P;
P1 = P(1:m,1:m)+1e-3;   % how to deal with 0 variance?
Hc = log2(2*pi*exp(1)*diag(P1))/2;

%-------------------- END --------------------

function sw = spectrum(w,obj,V)   % spectral density

fw = eval(exp(-1i*w),obj);
sw = fw*V*fw'/(2*pi);