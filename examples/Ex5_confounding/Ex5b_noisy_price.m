% This code solves the beauty-contest model in Section 2 of Huo and Pedroni
% (2021), assuming firms can observe a noisy past history of aggregate
% prices. 

clear
close all
clc 

% Parameters
alpha = 0.90;
rho = 0.90;
tau = 1;
kappa = 10; 

% Variables
xp_lag  = 1;    % lagged aggregate price: p_{t-1} 
xpi     = 2;    % individual optimal price: p_{i,t}
xzi     = 3;    % individual noisy signal on lagged price 

% Shocks
sq = 1;          % nominal aggregate demand shock 
ss = 2;          % idiosyncratic private signal about aggregate demand
sv = 3;          % idiosyncratic noise on lagged price 

% Innovations
ee = 1;
eu = 2;
ev = 3; 

% Dimensions
nx = 3;
ns = 3;
ne = 3; 

% Model
Ax0 = zeros(nx);
Aa0 = zeros(nx);  
Aa1 = zeros(nx);  
As0 = zeros(nx,ns); 
Ba0 = zeros(nx);  
Bs0 = zeros(nx,ns);
Cs1 = zeros(ns);
Ds0 = zeros(ns,ne);
Ds1 = zeros(ns,ne);
V = zeros(ne); 

% Eqn 1: Firm's pricing decision 
Ax0(1,xpi) = 1;
Ba0(1,xpi) = -alpha;
Bs0(1,sq) = -(1-alpha);

% Eqn 2: Define lagged price 
Ax0(2,xp_lag) = 1;
Aa1(2,xpi) = -1; 

% Eqn 3: Define noisy price signal 
Ax0(3,xzi) = 1;
Ax0(3,xp_lag) = -1; 
As0(3,sv) = -1; 

% Specify shocks 
% AR(1) aggregate demand shock 
Cs1(sq,sq) = rho;
Ds0(sq,ee) = 1;

% private signal structure
Cs1(ss,ss) = rho;
Ds0(ss,ee) = 1;
Ds0(ss,eu) = 1;
Ds1(ss,eu) = -rho; 

% noise on lagged price 
Ds0(sv,ev) = 1; 

V(ee,ee) = 1;
V(eu,eu) = 1/tau; 
V(ev,ev) = 1/kappa; 

agg = {xp_lag,ee};          % {var_id, inn_id}
sig = {
    1, xp_lag, ss, false    % {eqn_id, end_id, exo_id, ave_ex}
    };

% Solve model
% Solve model
m = ztran(nx,ns);
m.Ax = {Ax0};
m.Aa = {Aa0, Aa1};
m.As = {As0};
m.Ba = {Ba0};
m.Bs = {Bs0};
m.C = {Cs1};
m.D = {Ds0,Ds1};
m.V = V;
m.agg = agg;
m.sig = sig; 

m = solve(m,'dft',1000,'crit',1e-8,'nit',{10,2000,10},'step',0.1,...
    'arma', {5,5,false});

% IRF
T = 40;
var_lab = {'$p_{t-1}$', '$p_{t,i}$', '$z_{t,i}$'};
inn_lab = {'$\eta_t$','$u_{t,i}$', '$\nu_{t,i}$'};

figure('Name','IRF')
for i = 1:ne
    imp = zeros(ne,T);
    imp(i,1) = 1;
    res = irf(m.sol,imp);
    if i==eu
        res([xp_lag xzi],:) = 0;
    elseif i==ev
        res(xp_lag,:) = 0;
    end
    for j = 1:nx
        subplot(ne,nx,(i-1)*nx+j)
        plot(1:T,res(j,:),'linewidth',1.2)
        xlim([1 T])
        xticks([1 round(T/2) T])
        if i==1
            title(var_lab{j},'Interpreter','latex','FontSize',12)
        end
        if j==1
            ylabel(inn_lab{i},'Interpreter','latex','FontSize',12)
        end
    end
end

