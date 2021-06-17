function [U,S,V] = svdc(A,crit)
% Function SVDC
%
% Purpose:    Singular value decomposition in 'c'ompact form
%
% Format:     [U,S,V] = svdc(A)
%
% Input:      A         matrix to be decomposed (A = U*S*V')
%             crit      criterion for nonzero eigenvalues
%
% Output:     U         U'*U = I, not necessarily square
%             S         square and diagonal
%             V         V'*V = I, not necessarily square
%
% Adapted from Chris Sims' gensys.m
% Updated: September 14, 2019

if nargin<2
    crit = 1e-6;
end
[U,S,V] = svd(A);
ms = min(size(S));
bigev = diag(S(1:ms,1:ms))>crit;
U = U(:,bigev);
V = V(:,bigev);
S = S(bigev,bigev);