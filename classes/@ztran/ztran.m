classdef ztran < handle
% Handle class ZTRAN
% Written by Fei Tan, Saint Louis University
% Updated: August 22, 2023

    %% -------------------------------------------
    %                 Properties
    %---------------------------------------------
    
    properties
        % System matrix coeffs (cell); current & lag
        Ax         % {Ax0, Ax1, ...} - endogenous variable (cell)
        Aa         % {Aa0, Aa1, ...} - aggregate x (cell)
        As         % {As0, As1, ...} - exogenous shock (cell)
        
        % System matrix coeffs (cell); current & future
        Bx         % {Bx0, Bx1, ...} - expected endogenous variable (cell)
        Ba         % {Ba0, Ba1, ...} - expected aggregate x (cell)
        Bs         % {Bs0, Bs1, ...} - expected exogenous shock (cell)
        
        % Exogenous shock VARMA representation
        C          % {C1, C2, ...} - VAR matrix coeffs (cell)
        D          % {D0, D1, ...} - VMA matrix coeffs (cell)
        V          % innovation cov matrix
        
        % Additional restrictions
        agg        % {var_id, inn_id} - aggregate variable & innovation (cell)
        sig        % {eqn_id, end_id, exo_id, ave_ex} - signal structure (cell)
    end
    
    properties (GetAccess = public, SetAccess = protected)
        % Model solution
        apf        % {z(k), fz(:,:,k)} - analytic policy function (cell)
        sol        % varma object
                   % sol.AR = {AR1, AR2, ...} - VAR matrix coeffs (cell)
                   % sol.MA = {MA0, MA1, ...} - VMA matrix coeffs (cell)
        retcode    % diagnostic return code
                   % 0 - normal exit
                   % 1 - no solution
                   % 2 - more iterations needed
    end
    
    %% -------------------------------------------
    %                   Methods
    %---------------------------------------------
    
    methods
        % Set up model
        function obj = ztran(nx,ns)
            obj.Aa = {zeros(nx)};      % nx = number of variables
            obj.As = {zeros(nx,ns)};   % ns = number of shocks
            obj.Bx = {zeros(nx)};
            obj.Ba = {zeros(nx)};
            obj.Bs = {zeros(nx,ns)};
            obj.C = {zeros(ns)};
        end
        
        function set.Ax(obj,Ax)
            if ~iscell(Ax) || isempty(Ax)
                error('Ax must be non-empty cell array.')
            end
            obj.Ax = Ax;
        end
        
        function set.Aa(obj,Aa)
            if ~iscell(Aa) || isempty(Aa)
                error('Aa must be non-empty cell array.')
            end
            obj.Aa = Aa;
        end
        
        function set.As(obj,As)
            if ~iscell(As) || isempty(As)
                error('As must be non-empty cell array.')
            end
            obj.As = As;
        end
        
        function set.Bx(obj,Bx)
            if ~iscell(Bx) || isempty(Bx)
                error('Bx must be non-empty cell array.')
            end
            obj.Bx = Bx;
        end
        
        function set.Ba(obj,Ba)
            if ~iscell(Ba) || isempty(Ba)
                error('Ba must be non-empty cell array.')
            end
            obj.Ba = Ba;
        end
        
        function set.Bs(obj,Bs)
            if ~iscell(Bs) || isempty(Bs)
                error('Bs must be non-empty cell array.')
            end
            obj.Bs = Bs;
        end
        
        function set.C(obj,C)
            if ~iscell(C) || isempty(C)
                error('C must be non-empty cell array.')
            end
            obj.C = C;
        end
        
        function set.D(obj,D)
            if ~iscell(D) || isempty(D)
                error('D must be non-empty cell array.')
            end
            obj.D = D;
        end
        
        function set.agg(obj,agg)
            if ~iscell(agg) || isempty(agg)
                error('agg must be non-empty cell array.')
            end
            obj.agg = agg;
        end
        
        function set.sig(obj,sig)
            if ~iscell(sig) || isempty(sig)
                error('sig must be non-empty cell array.')
            end
            obj.sig = sig;
        end
        
        % Solve model
        solve(obj,varargin)
    end
end