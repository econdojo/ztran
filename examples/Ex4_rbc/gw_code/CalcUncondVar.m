function [Z_sigma] = CalcUncondVar(A, B,Q)

% calculates the unconditional variance of a system
% Z(t) = A Z(t-1) + B v(t)
% where A is a matrix, v is a vector of shocks 

% the var / covar matrix is given by
% vec(sigma) = (I - A xx A) *[B xx B] * vec(Q) NB xx is Kronecker product
pack
A_kron = kron(A, A);
B_kron = kron(B, B);
VecQ = reshape(Q,size(Q,1)*size(Q,2),1);
% the "\" operator speeds up the calc of the inverse - not sure how
% in general this is horribly inefficient as Z can be big, so kron (z) is enormous
Temp =  (eye(size(A_kron))-A_kron)\B_kron*VecQ;
Z_sigma=reshape(Temp,size(A,1),size(A,2));




