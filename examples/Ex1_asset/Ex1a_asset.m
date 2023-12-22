% This code solves a simple asset pricing model that permits an analytical
% solution. See Section 2 of Tan & Wu (2021): Analytic Policy Function
% Iteration.

clear
close all
clc

% Parameters
beta = 0.98;
gamma = 1.5;
rho = 0.9;
sig_e = 1;
sig_v = 3;

% Variables
xp = 1;

% Shocks
sd = 1;
ss = 2;

% Innovations
ee = 1;
ev = 2;

% Dimensions
nx = 1;
ns = 2;
ne = 2;

% Model
Ax0 = zeros(nx);
As0 = zeros(nx,ns);
Bx0 = zeros(nx);
Bx1 = zeros(nx);
C1 = zeros(ns);
D0 = zeros(ns,ne);
D1 = zeros(ns,ne);
V = zeros(ne);

Ax0(1,xp) = 1;
As0(1,sd) = -1;
Bx1(1,xp) = -beta;

C1(sd,sd) = rho;
D0(sd,ee) = 1;
D1(sd,ee) = -gamma;
D0(ss,ee) = 1;
D0(ss,ev) = 1;

V(ee,ee) = sig_e^2;
V(ev,ev) = sig_v^2;

agg = {xp,ee};     % {var_id, inn_id}
sig = {
    1, xp, ss, true     % {eqn_id, end_id, exo_id, ave_ex}
    };

% Solve model
m = ztran(nx,ns);
m.Ax = {Ax0};
m.As = {As0};
m.Bx = {Bx0,Bx1};
m.C = {C1};
m.D = {D0,D1};
m.V = V;
m.agg = agg;
m.sig = sig;
m.solve('dft',1000,'crit',1e-8,'nit',{10,500,10});

% IRF
T = 40;
var_lab = {'$p$'};
inn_lab = {'$\varepsilon$','$\nu_{i}$'};

figure('Name','IRF')
for i = 1:ne
    imp = zeros(ne,T);
    imp(i,1) = 1;
    res = m.sol.irf(imp);
    if ismember(i,setdiff(1:ne,agg{2}))
        res(agg{1},:) = 0;
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
