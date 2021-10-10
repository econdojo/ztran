% This code solves a prototypical DSGE model in Section 5.1.1 of Han, Tan &
% Wu (2021): Analytic Policy Function Iteration.

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

rho_a = 0.9764;     % AR(1) of technology shock 
rho_pi = 0.3411;    % AR(1) of inflation measurement error shock 
rho_y = 0.9541;     % AR(1) of output measurement error shock 
rho_m = 0.9468;     % AR(1) of monetary policy shock 

sig_a = 1.4208;     % 100*std. dev. of aggregate technology innovations 
sig_ai = 2.6068;    % 100*std. dev. of idiosyncratic technology innovations
sig_pi = 0.2686;    % 100*std. dev. of inflation measurement error innovations
sig_y = 1.0448;     % 100*std. dev. of output measurement error innovations
sig_m = 0.8474;     % 100*std. dev. of monetary policy innovations 

% Variables
% aggregate endogenous variables 
xy  = 1; 
xpi = 2; 
xi  = 3;
% individual endogenous variables 
xzi = 4;

% Shocks 
sa  = 1;    % technology shock    
spi = 2;    % MP inflation measurement error shock  
sy  = 3;    % MP output measurement error shock
sm  = 4;    % MP monetary policy shock 
sai = 5; 

% Innovations
% aggregate innovations
ea  = 1;
epi = 2; 
ey  = 3;
em  = 4; 
% idiosyncratic innovations 
eai = 5;

% Dimensions
nx = 4;
ns = 5;
ne = 5; 

% Model
Ax0 = zeros(nx);
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

V(ea,  ea) = sig_a^2;
V(epi,epi) = sig_pi^2;
V(ey,  ey) = sig_y^2; 
V(em,  em) = sig_m^2;
V(eai,eai) = sig_ai^2;

agg = {[xy xpi xi],[ea epi ey em]};     % {var_id, inn_id}
sig = {                                 % {eqn_id, end_id, exo_id, ave_ex}
    1, [], [], false                  % full-information households 
    3, xi, sai, false                 % incomplete information firms 
    }; 

% Solve model
m = ztran(nx,ns);
m.Ax = {Ax0};
m.Aa = {Aa0};
m.As = {As0}; 
m.Bx = {Bx0,Bx1}; 
m.Bs = {Bs0};
m.C = {Cs1};
m.D = {Ds0,Ds1};
m.V = V;
m.agg = agg;
m.sig = sig; 

m = solve(m,'dft',1000,'crit',1e-8,'nit',{10,1000,10},'step',0.5,...
    'arma', {10,10,false});

% IRF
T = 20;
var_lab = {'$y$', '$\pi$', '$i$', '$z_i$'};
inn_lab = {'$\varepsilon_a$', '$\varepsilon_\pi$', '$\varepsilon_y$',...
    '$\varepsilon_m$', '$\nu_{i}$'};

figure('Name','IRF')
for i = 1:ne
    imp = zeros(ne,T);
    imp(i,1) = 1;
    res = irf(m.sol,imp);
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
