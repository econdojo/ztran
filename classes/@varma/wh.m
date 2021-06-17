function [fz1,fz2] = wh(z,h,obj1,obj2,V,N)
% Function WH
%
% Purpose:    Evaluate Wiener-Hopf prediction formula
%
% Format:     [fz1,fz2] = wh(z,h,obj1,obj2,V,N)
%
% Input:      z(k)       point of evaluation
%             h          h-step-ahead prediction (integer), E_t[y(t+h)]
%                        infinite sum discount factor (0<h<1), h^0 E_t[y(t)]+h^1 E_t[y(t+1)]+h^2 E_t[y(t+2)]+...
%             obj1       signal variable varma object
%             obj2       forecast variable varma object
%             V          innovation cov matrix
%             N          number of DFT points
%
% Output:     fz1(:,:,k) innovation MA evaluated at z(k)
%             fz2(:,:,k) signal MA evaluated at z(k)
%
% Written by Fei Tan, Saint Louis University
% Updated: September 30, 2020

%% -------------------------------------------
%           Wiener-Hopf Prediction
%---------------------------------------------

% Initialization
nz = length(z);
[m1,n1] = size(obj1.MA{1});
m2 = size(obj2.MA{1},1);
SSR = ss(obj1,V,1:m1);
ns = size(SSR.F,1);
fz1 = zeros(m2,n1,nz);
fz2 = zeros(m2,m1,nz);
N = N+mod(N,2);
gz = zeros(m2,m1,N);

% Evaluate annihilation operator
for k = 0:N-1
    zz = exp(-1i*2*pi*k/N);
    e1 = eval(zz,obj2);
    e2 = eval(zz,obj1);
    e3 = eye(m1)+SSR.H/(eye(ns)-SSR.F*zz)*SSR.K*zz;
    if h>0 && h<1
        gz(:,:,k+1) = e1*V*e2'/(e3')*(zz/(zz-h));% note: conjugate transpose
    else
        gz(:,:,k+1) = e1*V*e2'/(e3')*zz^(-h);
    end
end
MA = real(ifftshift(ifft(gz,[],3),3)); % inverse DFT
MA = MA(:,:,(N/2+1):end);              % remove negative powers
MA = squeeze(num2cell(MA,[1 2]));
obj3 = varma({},MA);

% Evaluate W-H prediction
for k = 1:nz
    e1 = eval(z(k),obj3);
    e2 = eye(m1)+SSR.H/(eye(ns)-SSR.F*z(k))*SSR.K*z(k);
    e3 = eval(z(k),obj1);
    [U,S,V] = svdc(SSR.V);
    fz2(:,:,k) = e1*(V/S*U')/e2;
    fz1(:,:,k) = fz2(:,:,k)*e3;
end

%-------------------- END --------------------
