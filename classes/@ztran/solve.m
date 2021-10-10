function obj = solve(obj,varargin)
% Function SOLVE
%
% Purpose:    Implement analytic policy function iteration (APFI) algorithm
%
% Format:     obj = solve(obj,varargin)
%
% Input:      obj       initial ztran object
%             varargin  (optional) string/value pair
%                       'apf'     - {z(k), fz(:,:,k)}
%                                   z(:) points of evaluation in (-1,1)
%                                   fz(:,:,:) initial function values
%                       'nit'     - {min_iter, max_iter, disp_iter} ({10, 100, 0})
%                       'crit'    - convergence criterion (1e-5)
%                       'dft'     - number of DFT points (5000)
%                       'arma'    - {AR order, MA order, reduction} ({5, 5, false})
%                       'step'    - step size (1)
%
% Output:     obj       updated ztran object
%
% Written by Fei Tan, Saint Louis University
% Updated: September 30, 2020

%% -------------------------------------------
%               APFI Algorithm
%---------------------------------------------

% Initialization
lx = length(obj.Ax);
nx = size(obj.Ax{1},1);
[ns,ne] = size(obj.D{1});
nex = size(obj.sig,1);  % number of expectational equations
ind_var = setdiff(1:nx,obj.agg{1});
ind_inn = setdiff(1:ne,obj.agg{2});
Ax = cell2mat(obj.Ax(:));
ind_eqn = any(reshape(Ax(:,ind_var),nx,[]),2);
obj.retcode = 2;

% Optional inputs
z = linspace(-0.99,0.99,50);
nz = length(z);
fz = zeros(nx,ne,nz);
nit1 = 10; nit2 = 1000; nit3 = 0;
crit = 1e-6;
N = 1000;
p = 5; q = 5; or = false;
step = 1;

narg = length(varargin);
for k = 1:2:narg
    switch varargin{k}
        case 'apf', z = varargin{k+1}{1}; fz = varargin{k+1}{2}; nz = length(z);
            if any(fz(obj.agg{1},ind_inn,:),'all')
                error('fz must be 0 for aggregate variables w.r.t. individual innovations.')
            end
        case 'nit', nit1 = varargin{k+1}{1}; nit2 = varargin{k+1}{2}; nit3 = varargin{k+1}{3};
        case 'crit', crit = varargin{k+1};
        case 'dft', N = varargin{k+1}; % N ~ 2*log(1e-5)/log(max_root)
        case 'arma', p = varargin{k+1}{1}; q = varargin{k+1}{2}; or = varargin{k+1}{3};
        case 'step', step = varargin{k+1};
        otherwise
            error('Unrecognized optional argument.')
    end
end

% Lagged (x,a,s)
z = reshape(z,1,1,nz);
la = length(obj.Aa);
lxa = max(lx,la);
Lx = zeros(nx,ne,nz,lxa);
Lx(:,:,:,1) = fz;
ls = length(obj.As);
Ls = zeros(ns,ne,nz,ls);
Ls(:,:,:,1) = eval(z,varma(obj.C,obj.D));
if ls>1
    zz = repmat(z,ns,ne);
    for k = 2:ls
        Ls(:,:,:,k) = Ls(:,:,:,k-1).*zz;
    end
end
if lxa>1
    zz = repmat(z,nx,ne);
end

% Expected (x,a,s)
hx = length(obj.Bx);
ha = length(obj.Ba);
hs = length(obj.Bs);
Ex = zeros(nx,ne,nz,hx,nex);
Ea = zeros(nx,ne,nz,ha,nex);
Es = zeros(ns,ne,nz,hs,nex);

% APFI
for iter = 1:nit2
    % Lagged (x,a)
    for k = 2:lxa
        Lx(:,:,:,k) = Lx(:,:,:,k-1).*zz;
    end
    
    % Expected (x,a,s)
    for j = 1:nex
        sigx = varma.fit(z,Lx(obj.sig{j,2},:,:,1),p,q,or); % endogenous signal
        sigs = varma.fit(z,Ls(obj.sig{j,3},:,:,1),p,q,or); % exogenous signal
        sig = sigx+sigs;     % joint signal
        if isempty(sig.MA)   % complete info
            sig = varma({zeros(ne)},{eye(ne)});
        end
        
        for k = 1:hx
            if any(obj.Bx{k}(obj.sig{j,1},:))
                predx = varma.fit(z,Lx(obj.Bx{k}(obj.sig{j,1},:)~=0,:,:,1),p,q,or);
                [Ex(obj.Bx{k}(obj.sig{j,1},:)~=0,:,:,k,j),~] = wh(z,k-1,sig,predx,obj.V,N);
                if obj.sig{j,4}   % average expectation
                    Ex(obj.Bx{k}(obj.sig{j,1},:)~=0,ind_inn,:,k,j) = 0;
                end
            end
        end
        
        for k = 1:ha
            if any(obj.Ba{k}(obj.sig{j,1},:))
                La = Lx(obj.Ba{k}(obj.sig{j,1},:)~=0,:,:,1);
                La(:,ind_inn,:) = 0;
                preda = varma.fit(z,La,p,q,or);
                [Ea(obj.Ba{k}(obj.sig{j,1},:)~=0,:,:,k,j),~] = wh(z,k-1,sig,preda,obj.V,N);
                if obj.sig{j,4}
                    Ea(obj.Ba{k}(obj.sig{j,1},:)~=0,ind_inn,:,k,j) = 0;
                end
            end
        end
        
        for k = 1:hs
            if any(obj.Bs{k}(obj.sig{j,1},:))
                preds = varma.fit(z,Ls(obj.Bs{k}(obj.sig{j,1},:)~=0,:,:,1),p,q,or);
                [Es(obj.Bs{k}(obj.sig{j,1},:)~=0,:,:,k,j),~] = wh(z,k-1,sig,preds,obj.V,N);
                if obj.sig{j,4}
                    Es(obj.Bs{k}(obj.sig{j,1},:)~=0,ind_inn,:,k,j) = 0;
                end
            end
        end
    end
    
    % Update policy function
    res = 0;
    for k = 1:nz
        tmp = zeros(nx,ne);
        for j = 1:la
            La = Lx(:,:,k,j);
            La(:,ind_inn) = 0;
            tmp = tmp+obj.Aa{j}*La;
        end
        
        for j = 1:ls
            tmp = tmp+obj.As{j}*Ls(:,:,k,j);
        end
        
        for j = 1:hx
            for i = 1:nex
                tmp(obj.sig{i,1},:) = tmp(obj.sig{i,1},:)+obj.Bx{j}(obj.sig{i,1},:)*Ex(:,:,k,j,i);
            end
        end
        
        for j = 1:ha
            for i = 1:nex
                tmp(obj.sig{i,1},:) = tmp(obj.sig{i,1},:)+obj.Ba{j}(obj.sig{i,1},:)*Ea(:,:,k,j,i);
            end
        end
        
        for j = 1:hs
            for i = 1:nex
                tmp(obj.sig{i,1},:) = tmp(obj.sig{i,1},:)+obj.Bs{j}(obj.sig{i,1},:)*Es(:,:,k,j,i);
            end
        end
        
        Ax = obj.Ax{1};
        for j = 2:lx
            Ax = Ax+obj.Ax{j}*z(k)^(j-1);
        end
        
        fz(ind_var,ind_inn,k) = -Ax(ind_eqn,ind_var)\tmp(ind_eqn,ind_inn);
        fz(:,obj.agg{2},k) = -Ax\tmp(:,obj.agg{2});
        res = max(res,max(max(abs(Ax*fz(:,:,k)+tmp))));
    end
    
    % Check convergence
    fz = step*fz+(1-step)*Lx(:,:,:,1);
    gap = max(max(max(abs(fz-Lx(:,:,:,1)))));
    gap = gap/(max(max(max(abs(Lx(:,:,:,1)))))+1e-5);
    if ~mod(iter,nit3)
        disp(['Iter ' num2str(iter) ': gap = ' num2str(gap) ', res = ' num2str(res)]); disp(' ');
    end
    
    if gap<crit && iter>=nit1% warm-up
        if res<1e-5
            obj.retcode = 0;
        else
            obj.retcode = 1;
        end
        break
    else
        Lx(:,:,:,1) = fz;
    end
end

% Output
obj.apf = {z(:),fz};
obj.sol = varma.fit(z,fz,p,q,or);

if nit3
    switch obj.retcode
        case 0
            fprintf('Iter %d: APFI converged, solution exists.\n\n',iter);
        case 1
            fprintf('Iter %d: APFI converged, NO solution exists.\n\n',iter);
        case 2
            fprintf('Iter %d: APFI NOT converged.\n\n',nit2);
    end
end

%-------------------- END --------------------
