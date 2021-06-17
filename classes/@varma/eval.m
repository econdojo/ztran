function fz = eval(z,obj)
% Function EVAL
%
% Purpose:    Evaluate z-transform of VARMA representation
%
% Format:     fz = eval(z,obj)
%
% Input:      z(k)      point of evaluation
%             obj       varma object
%
% Output:     fz(:,:,k) VARMA representation evaluated at z(k)
%
% Written by Fei Tan, Saint Louis University
% Updated: September 30, 2020

%% -------------------------------------------
%              VARMA Evaluation
%---------------------------------------------

% Initialization
nz = length(z);
p = length(obj.AR);
q = length(obj.MA)-1;
[m,n] = size(obj.MA{1});
fz = zeros(m,n,nz);

% Evaluation
for k = 1:nz
    AA = eye(m);
    BB = obj.MA{1};
    for j = 1:p
        AA = AA-obj.AR{j}*z(k)^j;
    end
    for j = 1:q
        BB = BB+obj.MA{j+1}*z(k)^j;
    end
    fz(:,:,k) = AA\BB;
end

%-------------------- END --------------------
