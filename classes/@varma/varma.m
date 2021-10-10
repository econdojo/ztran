classdef varma
% Class VARMA
% Written by Fei Tan, Saint Louis University
% Updated: September 30, 2020

    %% -------------------------------------------
    %                 Properties
    %---------------------------------------------

    properties
        % VARMA process: x(t) = AR1*x(t-1) + ... + ARp*x(t-p) + MA0*e(t) + ... + MAq*e(t-q)
        AR    % {AR1, AR2, ... ARp} - VAR matrix coeffs (cell)
        MA    % {MA0, MA1, ... MAq} - VMA matrix coeffs (cell)
    end
    
    %% -------------------------------------------
    %                   Methods
    %---------------------------------------------
    
    methods
        % Constructor
        function obj = varma(AR,MA)
            if ~iscell(AR) || ~iscell(MA)
                error('AR/MA must be cell array.')
            end
            obj.AR = AR;
            obj.MA = MA;
        end
        
        % Combine two VARMAs with same innovations
        obj = plus(obj1,obj2)
        
        % Evaluate z-transform of VARMA representation
        fz = eval(z,obj)
        
        % Represent VARMA process in state space form
        SSR = ss(obj,V,sid)
        
        % Evaluate Wiener-Hopf prediction formula
        [fz1,fz2] = wh(z,h,obj1,obj2,V,N)
        
        % Simulate VARMA impulse response function
        res = irf(obj,imp,crit)
        
        % Compute Shannon information flow
        [Hu,Hc,P0,P1,P2] = info(obj,V,sid)
    end
    
    methods (Static)
        % Fit VARMA representation to data points
        obj = fit(z,fz,p,q,or)
    end
end