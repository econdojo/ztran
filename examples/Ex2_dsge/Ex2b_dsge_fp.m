% This code solves a DSGE model with both rule-based monetary and fiscal
% policies. See Section 5.1.2 of Tan & Wu (2021): Analytic Policy
% Function Iteration.

clear
close all
clc

% Parameters
beta = 0.99;        % quarterly discount rate
theta = 0.3608;     % price stickiness
gamma = 1;          % relative risk aversion
kappa = 0;          % inverse of Frisch elasticity

phi_pi = 1.6782;    % MP response to inflation
phi_y = 0.6731;     % MP response to output

phi_b = 1.50;       % FP response to lagged debt

rho_a = 0.9764;     % AR(1) of technology shock
rho_pi = 0.3411;    % AR(1) of inflation measurement error shock
rho_y = 0.9541;     % AR(1) of output measurement error shock
rho_m = 0.9468;     % AR(1) of monetary policy shock
rho_f = 0.51;       % AR(1) of fiscal policy shock

sig_a = 1.4208;     % 100*std. dev. of aggregate technology innovations
sig_ai = 2.6068;    % 100*std. dev. of idiosyncratic technology innovations
sig_pi = 0.2686;    % 100*std. dev. of inflation measurement error innovations
sig_y = 1.0448;     % 100*std. dev. of output measurement error innovations
sig_m = 0.8474;     % 100*std. dev. of monetary policy innovations
sig_f = 0.43;       % 100*the std of FP innovation

% Variables
% aggregate endogenous variables
xy  = 1;
xpi = 2;
xi  = 3;
xs  = 4;    % primary surplus
xb  = 5;    % debt
% individual endogenous variables
xzi = 6;

% Shocks
sa  = 1;    % technology shock
spi = 2;    % MP inflation measurement error shock
sy  = 3;    % MP output measurement error shock
sm  = 4;    % MP monetary policy shock
sf  = 5;    % FP fiscal policy shock
sai = 6;

% Innovations
% aggregate innovations
ea  = 1;
epi = 2;
ey  = 3;
em  = 4;
ef  = 5;
% idiosyncratic innovations
eai = 6;

% Dimensions
nx = 6;
ns = 6;
ne = 6;

% Model
Ax0 = zeros(nx);
Ax1 = zeros(nx);
Aa0 = zeros(nx);
As0 = zeros(nx,ns);
Bx0 = zeros(nx);
Bx1 = zeros(nx);
Bs0 = zeros(nx,ns);
Cs1 = zeros(ns);
Ds0 = zeros(ns,ne);
Ds1 = zeros(ns,ne);
V = zeros(ne);

% Eqn 1: IS equation
Ax0(1,xy) = 1;
Ax0(1,xi) = 1/gamma;
Bx1(1,xy) = -1;
Bx1(1,xpi) = -1/gamma;

% Eqn 2: monetary policy rule
Ax0(2,xi) = 1;
Ax0(2,xpi) = -phi_pi;
Ax0(2,xy) = -phi_y;
As0(2,spi) = -phi_pi;
As0(2,sy) = -phi_y;
As0(2,sm) = -1;

% Eqn 3: Phillips Curve
Ax0(3,xzi) = 1;
As0(3,sai) = (1-beta*theta);
Bs0(3,sa) = (1-beta*theta)*kappa;
Bx0(3,xy) = -(kappa+gamma)*(1-beta*theta);
Bx0(3,xpi) = -1;
Bx1(3,xzi) = -beta*theta;

% Eqn 4: price aggregation
Ax0(4,xpi) = 1;
Aa0(4,xzi) = -(1-theta);

% Eqn 5: Fiscal policy
Ax0(5,xs) = 1;
Ax1(5,xb) = -phi_b;
As0(5,sf) = -1;

% Eqn 6: Government budget constraint
Ax0(6,xb) = 1;
Ax1(6,xi) = -1/beta;
Ax1(6,xb) = -1/beta;
Ax0(6,xpi) = 1/beta;
Ax0(6,xs) = 1/beta-1;

% Specify shocks
% private signal structure
Cs1(sai,sai) = rho_a;
Ds0(sai,ea) = 1;
Ds0(sai,eai) = 1;
Ds1(sai,eai) = -rho_a;

% AR(1) aggregate technology shock
Cs1(sa,sa) = rho_a;
Ds0(sa,ea) = 1;

% AR(1) inflation measurement error shock
Cs1(spi,spi) = rho_pi;
Ds0(spi,epi) = 1;

% AR(1) output measurement error shock
Cs1(sy,sy) = rho_y;
Ds0(sy,ey) = 1;

% AR(1) monetary policy shock
Cs1(sm,sm) = rho_m;
Ds0(sm,em) = 1;

% AR(1) fiscal policy shock
Cs1(sf,sf) = rho_f;
Ds0(sf,ef) = 1;

V(ea,  ea) = sig_a^2;
V(epi,epi) = sig_pi^2;
V(ey,  ey) = sig_y^2;
V(em,  em) = sig_m^2;
V(ef,  ef) = sig_f^2;
V(eai,eai) = sig_ai^2;

agg = {[xy xpi xi xs xb],[ea epi ey em ef]};  % {var_id, inn_id}

% Consider four cases of different incomplete information sets.
sig_case = 1;

if sig_case == 1
    % Case 1
    sig = {                              % {eqn_id, end_id, exo_id, ave_ex}
        1, xi, [], false
        3, [xi, xs], sai, false
        };
elseif sig_case == 2
    % Case 2
    sig = {
        1, xi, [], false
        3, [xi, xb], sai, false
        };
elseif sig_case == 3
    % Case 3
    sig = {
        1, [xi, xs], [], false
        3, xi, sai, false
        };
elseif sig_case == 4
    % Case 4
    sig = {
        1, [xi, xb], [], false
        3, xi, sai, false
        };
end

% Solve model
m = ztran(nx,ns);
m.Ax = {Ax0, Ax1};
m.Aa = {Aa0};
m.As = {As0};
m.Bx = {Bx0,Bx1};
m.Bs = {Bs0};
m.C = {Cs1};
m.D = {Ds0,Ds1};
m.V = V;
m.agg = agg;
m.sig = sig;

m.solve('dft',1000,'crit',1e-8,'nit',{10,2000,10},'step',0.1,...
    'arma', {11,11,false});

% IRF
T = 20;
var_lab = {'$y$', '$\pi$', '$i$', '$s$', '$b$', '$z_i$'};
inn_lab = {'$\varepsilon_a$', '$\varepsilon_\pi$', '$\varepsilon_y$',...
    '$\varepsilon_m$', '$\varepsilon_f$', '$\nu_{i}$'};

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
