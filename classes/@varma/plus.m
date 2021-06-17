function obj = plus(obj1,obj2)
% Function PLUS
%
% Purpose:    Combine two VARMAs with same innovations
%
% Format:     [obj] = plus(obj1,obj2)
%
% Input:      obj1      varma object 1
%             obj2      varma object 2
%
% Output:     obj       joint varma object
%
% Written by Fei Tan, Saint Louis University
% Updated: September 30, 2020

%% -------------------------------------------
%                 Joint VARMA
%---------------------------------------------

% Initialization
if isempty(obj1.MA)
    obj = obj2;
    return
elseif isempty(obj2.MA)
    obj = obj1;
    return
end
p1 = length(obj1.AR);
q1 = length(obj1.MA);
[nx1,ne1] = size(obj1.MA{1});
p2 = length(obj2.AR);
q2 = length(obj2.MA);
[nx2,ne2] = size(obj2.MA{1});

% Joint VAR
r = max(p1,p2);
AR = cell(r,1);
for k = 1:r
    if k<=p1
        AR1 = obj1.AR{k};
    else
        AR1 = zeros(nx1);
    end
    if k<=p2
        AR2 = obj2.AR{k};
    else
        AR2 = zeros(nx2);
    end
    AR{k} = [AR1 zeros(nx1,nx2);zeros(nx2,nx1) AR2];
end

% Joint VMA
r = max(q1,q2);
MA = cell(r,1);
for k = 1:r
    if k<=q1
        MA1 = obj1.MA{k};
    else
        MA1 = zeros(nx1,ne1);
    end
    if k<=q2
        MA2 = obj2.MA{k};
    else
        MA2 = zeros(nx2,ne2);
    end
    MA{k} = [MA1;MA2];
end

% Joint VARMA
obj = varma(AR,MA);

%-------------------- END --------------------
