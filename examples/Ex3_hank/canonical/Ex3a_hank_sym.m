% This code solves a HANK model with symmetric incomplete information in
% Section 5.2 of Han, Tan & Wu (2021): Analytic Policy Function Iteration.

clear
close all
clc

% Parameters
beta = 0.99;        % quarterly discount rate
m1 = 0.55;          % MPC of group 1
m2 = 0.05;          % MPC of group 2
phi1 = 1.99; 2;     % group 1's exposure to business cycle
phi2 = 0.01; 0;     % group 2's exposure to business cycle
pi1 = 0.5;          % group 1's mass
pi2 = 0.5;          % group 2's mass

w1 = (1-m1)/beta;   % group 1's survival rate each period
w2 = (1-m2)/beta;   % group 2's survival rate each period

rho_r = 0.9;        % AR(1) persistence of real interest rate
sig_r = 1;          % std. dev. of real interest rate shock
sig_eta = 1.3801;   % std. dev. of private signal shock

sig_xi = 1;         % std. dev. of private endogenous signal

% Variables aggregate endogenous variables
xy  = 1;    % aggregate real output
xc1 = 2;    % group 1 aggregate consumption
xc2 = 3;    % group 2 aggregate consumption
xs1 = 4;    % group 1 aggregate saving
xs2 = 5;    % group 2 aggregate saving

% individual endogenous variables
xx1i = 6;   % group 1's dummy to present value of real interest rate
xx2i = 7;   % group 2's dummy to present value of real interest rate
xz1i = 8;   % group 1's dummy to present value of real output
xz2i = 9;   % group 2's dummy to present value of real output

% Shocks
sr  = 1;    % shock to real interest rate
s1i = 2;    % shock to group 1's individual private signal
s2i = 3;    % shock to group 2's individual private signal

% Innovations aggregate innovations
er  = 1;

% idiosyncratic innovations
e1i = 2;
e2i = 3;

% Dimensions
nx = 9;
ns = 3;
ne = 3;

% Model
Ax0 = zeros(nx);
Ax1 = zeros(nx);
Aa0 = zeros(nx);
As0 = zeros(nx,ns);
As1 = zeros(nx,ns);
Bx0 = zeros(nx);
Bx1 = zeros(nx);
Bs0 = zeros(nx,ns);
Cs1 = zeros(ns);
Ds0 = zeros(ns,ne);
Ds1 = zeros(ns,ne);
V = zeros(ne);

% Eqn 1: group 1's aggregate consumption
Ax0(1, xc1) = 1;
Ax1(1, xs1) = -(1-beta*w1)/beta;
Aa0(1, xx1i) = 1;
Aa0(1, xz1i) = -1;

% Eqn 2: group 2's aggregate consumption
Ax0(2, xc2) = 1;
Ax1(2, xs2) = -(1-beta*w2)/beta;
Aa0(2, xx2i) = 1;
Aa0(2, xz2i) = -1;

% Eqn 3: group 1's budget constraint
Ax0(3, xc1) = 1;
Ax0(3, xs1) = 1;
Ax1(3, xs1) = -1/beta;
Ax0(3, xy) = -phi1;

% Eqn 4: group 2's budget constraint
Ax0(4, xc2) = 1;
Ax0(4, xs2) = 1;
Ax1(4, xs2) = -1/beta;
Ax0(4, xy) = -phi2;

% Eqn 5: market clearing condition
Ax0(5, xs1) = pi1;
Ax0(5, xs2) = pi2;

% Eqn 6: define dummy xx1i
Ax0(6, xx1i) = 1;
Bs0(6, sr) = -beta*w1;
Bx1(6, xx1i) = -beta*w1;

% Eqn 7: define dummy xx2i
Ax0(7, xx2i) = 1;
Bs0(7, sr) = -beta*w2;
Bx1(7, xx2i) = -beta*w2;

% Eqn 8: define dummy xz1i
Ax0(8, xz1i) = 1;
Bx0(8, xy) = -(1-beta*w1)*phi1;
Bx1(8, xz1i) = -beta*w1;

% Eqn 9: define dummy xz2i
Ax0(9, xz2i) = 1;
Bx0(9, xy) = -(1-beta*w2)*phi2;
Bx1(9, xz2i) = -beta*w2;

% Specify shocks AR(1) real interest rate shock
Cs1(sr,sr) = rho_r;
Ds0(sr,er) = 1;

% private signal structure
Cs1(s1i,s1i) = rho_r;
Ds0(s1i,er) = 1;
Ds0(s1i,e1i) = 1;
Ds1(s1i,e1i) = -rho_r;

Cs1(s2i,s2i) = rho_r;
Ds0(s2i,er) = 1;
Ds0(s2i,e2i) = 1;
Ds1(s2i,e2i) = -rho_r;


V(er,  er)  = sig_r^2;
V(e1i, e1i) = sig_eta^2;
V(e2i, e2i) = sig_eta^2;

agg = {[xc1 xc2 xs1 xs2 xy],er};   % {var_id, inn_id}

sig = {                              % {eqn_id, end_id, exo_id, ave_ex}
    6, [], s1i, false
    7, [], s2i, false
    8, [], s1i, false
    9, [], s2i, false
    };

% Solve model
m = ztran(nx,ns);
m.Ax = {Ax0, Ax1};
m.Aa = {Aa0};
m.As = {As0, As1};
m.Bx = {Bx0,Bx1};
m.Bs = {Bs0};
m.C = {Cs1};
m.D = {Ds0,Ds1};
m.V = V;
m.agg = agg;
m.sig = sig;

m.solve('dft',500,'crit',1e-8,'nit',{10,2000,10},'step',1,...
    'arma', {3,3,false});

% IRFs
T = 40;
var_lab = {'$c_1$', '$c_2$', '$s_1$', '$s_2$', '$y$', ...
    '$x_1^i$', '$x_2^i$', '$z_1^i$', '$z_2^i$', '$m_1^i$'};
inn_lab = {'$\varepsilon_r$', '$\eta_1^i$', '$\eta_2^i$'};

figure('Name','IRF')
for i = 1:ne
    imp = zeros(ne,T);
    imp(i,1) = -1;
    res = m.sol.irf(imp);
    res(abs(res)<1e-5) = 0; 
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