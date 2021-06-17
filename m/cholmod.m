function [R, indef, E] = cholmod(A)
% CHOLMOD Modified Cholesky factorization
%  R = cholmod(A) returns the upper Cholesky factor of A (same as CHOL)
%  if A is (sufficiently) positive definite, and otherwise returns a 
%  modified factor R with diagonal enries >= sqrt(delta) and
%  offdiagonal entries <= beta in absolute value,
%  where delta and beta are defined in terms of size of diagonal and
%  offdiagonal entries of A and the machine precision; see below.
%  The idea is to ensure that E = A - R'*R is reasonably small if A is 
%  not too far from being indefinite.  If A is sparse, so is R.
%  The output parameter indef is set to 0 if A is sufficiently positive
%  definite and to 1 if the factorization is modified.
%
%  The point of modified Cholesky is to avoid computing eigenvalues, 
%  but for dense matrices, MODCHOL takes longer than calling the built-in 
%  EIG, because of the cost of interpreting the code, even though it
%  only has one loop and uses vector operations.  

%  reference: Nocedal and Wright, Algorithm 3.4 and subsequent discussion
%  (not Algorithm 3.5, which is more complicated)
%  original algorithm is due to Gill and Murray, 1974
%  written by M. Overton, overton@cs.nyu.edu, last modified Feb 2005

%  convenient to work with A = LDL' where D is diagonal, L is unit
%  lower triangular, and so R = (LD^(1/2))'

% if sum(sum(abs(A-A'))) > 0
%     error('A is not symmetric')
% end
% Line 27-29 modified (by Fei Tan, 12/7/2016) as below:
A = (A+A')/2; % trick for symmetry

% set parameters governing bounds on L and D (eps is machine epsilon)

n = length(A);
diagA = diag(A);
gamma = max(abs(diagA));             % max diagonal entry
xi = max(max(abs(A - diag(diagA))));  % max offidagonal entry
delta = eps*(max([gamma+xi, 1]));
beta = sqrt(max([gamma, xi/n, eps]));
indef = 0;

% initialize d and L

d = zeros(n,1);
%if issparse(A)
%    L = speye(n);  % sparse identity
%else
    L = eye(n);    
%end

% there are no inner for loops, everything implemented with
% vector operations for a reasonable level of efficiency

for j = 1:n
    K = 1:j-1;  % column index: all columns to left of diagonal
                % d(K) doesn't work in case K is empty
    djtemp = A(j,j) - L(j,K)*(d(K,1).*L(j,K)');   % C(j,j) in book
    if j < n
        I = j+1:n;  % row index: all rows below diagonal
        Ccol = A(I,j) - L(I,K)*(d(K,1).*L(j,K)');  % C(I,j) in book
        theta = max(abs(Ccol));
        % guarantees d(j) not too small and L(I,j) not too big
        % in sufficiently positive definite case, d(j) = djtemp
        d(j) = max([abs(djtemp), (theta/beta)^2, delta]);
        L(I,j) = Ccol/d(j);
    else
        d(j) = max([abs(djtemp), delta]);
    end
    if d(j) > djtemp  % A was not sufficiently positive definite
        indef = 1;
    end
end
% convert to usual output format: replace L by L*sqrt(D) and transpose
for j=1:n
    L(:,j) = L(:,j)*sqrt(d(j));   % L = L*diag(sqrt(d)) bad in sparse case
end;
R = L';
if nargout == 3
    E = A - R'*R;
end