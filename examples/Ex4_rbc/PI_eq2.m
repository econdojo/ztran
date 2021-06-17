% This code solves the model of Graham and Wright (2010, JME)  
% Incomplete Information Version Using closed form solution
% Solve for PI equilibrium 2 where individual consumption is stationary, no
% random walk, idiosyncratic capital has explosive roots at z= beta

clear
close all
clc

% Parameters
beta = 0.99;        % quarterly discount rate
gamma = 5;          % inverse of Frisch elasticity
alpha = 0.667;      % labor share of output 
delta = 0.025;      % depreciation rate of capital 

rho_a = 0.9;        % persistence of aggregate technology 
rho_z = 0.9;        % persistence of labor productivity 
sig_a = 0.7;        % std. dev. of aggregate technology shock
sig_z = 4.9;        % std. dev. of idiosyncratic labor produtivity shock

% Steady state 
H_ss = 1/3;         % steady-state labor supply 
K_ss = ((1-alpha)/(1/beta-(1-delta)))^(1/alpha)*H_ss; 
I_ss = delta*K_ss; 
Y_ss = K_ss^(1-alpha)*(H_ss^alpha); 
C_ss = Y_ss - I_ss; 

% composite parameters in the two equation cannonical form
a_w = (1/beta)*(alpha/(1-alpha))*(1-beta*(1-delta));
a_c = C_ss/K_ss;
a_n = (1-H_ss)/(gamma*H_ss);
g_1 = 1+a_n*(1-alpha*beta*(1-delta));
g_2 = beta*(1+(1-alpha)*a_n);
g_3 = alpha*(1+a_n)+alpha*beta*a_w*(1+a_n);
g_4 = alpha*beta*(1+a_n);
g_5 = alpha*a_n+alpha*beta*(a_c+a_w*a_n);
g_6 = alpha*beta*a_n;

%% Analytical solution, notation follows my notes

% Exogenous polynomials
Qz = [1+a_n*(1-alpha*beta*(1-delta))  -beta*(1+(1-alpha)*a_n)];
Pcz = [alpha*a_n+alpha*beta*(a_c+a_w*a_n)  -alpha*beta*a_n];
Paz = [alpha*(1+a_n)+alpha*beta*a_w*(1+a_n)  -alpha*beta*(1+a_n)];

Phiz_1 = conv([1 -1],Qz);
Phiz_2 = -(1-beta*(1-delta))*Pcz;

% Denominator polynomials of c_t and r_kt
Phiz = sum_poly_coeff(Phiz_1, Phiz_2);
r = roots(Phiz);

r_1 = r(1); % outside root
r_2 = r(2); % inside root

% leading coefficients of Phi_z in the denominator
Phi_1 = 1+a_n*(1-alpha*beta*(1-delta));

%% Solve for endogenous confounding roots lambda

% Wold-Related Matrix of Polynomials

% Assuming rho_a = rho_z = rho.
tau_a = sig_a^(2)/(sig_a^2+sig_z^2);
tau_b = (sig_a*sig_z)/(sig_a^2+sig_z^2);

F_121l = tau_b*[-1 0 1];
F_221l = conv([-1 1], [1-tau_a -tau_a]);
F_11rl_1 = tau_a*conv([1 0],[-r_2 1]);
F_11rl_2 = -(1-tau_a)*[-1 r_2];
F_11rl = sum_poly_coeff(F_11rl_1, F_11rl_2 );
F_12rl_1 = tau_b*conv([1 0],[-r_2 1]);
F_12rl_2 = tau_b*[-1 r_2];
F_12rl = sum_poly_coeff(F_12rl_1, F_12rl_2 );

Ml_1 = conv(F_221l,F_11rl);
Ml_2 = conv(F_12rl,F_121l);
Ml = sum_poly_coeff(Ml_1, -Ml_2);

h_r2 = (r_2-1)*polyval(Paz,r_2)/polyval(Pcz,r_2);
a_r2 = sig_a/(1-rho_a*r_2);
c_t = h_r2*a_r2;

fl_1 =sig_a*conv(Paz,[1 -1]);
fl_1 = conv(fl_1,[1-(1-r_2)*tau_a  -(r_2+(1-r_2)*tau_a)]);
fl_2 = c_t*tau_a*conv([-rho_a 1],Pcz);
fl_2 = conv(fl_2,[-r_2 1]);
fl = sum_poly_coeff(fl_1,fl_2);

R_l = roots(fl); % One of the inside roots at z=r_2
lambda = R_l(end); % endogenous non-invertible roots generating the confounding dynamics

%% Given lambda, define functional constant theta_1 and theta_2, and Wold matrix polynomial function F(z;lambda)
% Numerator polynomials of F(z;lambda)
F_11z = sum_poly_coeff(tau_a*lambda*[-lambda 1], -(1-tau_a)*[1 -lambda]);
F_12z = tau_b*sum_poly_coeff(lambda*[-lambda 1],[1 -lambda]);
F_22z = sum_poly_coeff((1-tau_a)*lambda*[-lambda 1], -(tau_a)*[1 -lambda]); 

% Functional constant
M_l = polyval(F_22z,1)*polyval(F_11z,r_2)-polyval(F_12z,r_2)*polyval(F_12z,1); % function of lambda
%M_l = polyval(Ml,lambda);
theta_1 = (c_t*lambda*(1-lambda*r_2)*polyval(F_22z,1))/M_l;
theta_2 = (c_t*lambda*(1-lambda*r_2)*polyval(F_12z,1))/M_l;

% Numerator Polynomial function of C(z) and Rk(z)
N_1z = sum_poly_coeff(theta_1*F_11z,-theta_2*F_12z);
N_2z = sum_poly_coeff(theta_1*F_12z,-theta_2*F_22z);

Num_caz_1 = conv(conv(N_1z,Qz),[-rho_a 1]);
Num_caz_2 = (1-beta*(1-delta))*conv(sig_a*Paz,lambda*[-lambda 1]);
Num_caz = sum_poly_coeff(Num_caz_1,-Num_caz_2);

r_num_caz = roots(Num_caz);
r_3 = r_num_caz(1);
r_4 = r_num_caz(3); % two remaining roots
C_f = Num_caz(1); % leading coefficients

Num_rkz_1 = lambda*sig_a*conv(conv([1 -1],Paz),[-lambda 1]);
Num_rkz_2 = conv(conv(Pcz,N_1z),[-rho_a 1]);
Num_rkz = sum_poly_coeff(Num_rkz_1,-Num_rkz_2);

r_num_rkz = roots(Num_rkz);
r_5 = r_num_rkz(1); % outside roots, one inside roots at z = lambda, one inside root at z = r_2
R_f = Num_rkz(1);

%% Function Value based on analytical solution
z = linspace(-0.99,0.99,100);
nz = length(z);
fz = zeros(2,1,nz); % function values using closed-form solution, aggregate consumption and capital return

for i = 1:nz
    fz(1,1,i) = (1/sig_a)*C_f*(z(i)-r_3)*(z(i)-r_4)/(Phi_1*lambda*(1-lambda*z(i))*(1-rho_a*z(i))*(z(i)-r_1));
    fz(2,1,i) = (1/sig_a)*R_f*(z(i)-r_5)*(z(i)-lambda)/(Phi_1*lambda*(1-lambda*z(i))*(1-rho_a*z(i))*(z(i)-r_1));
end

ccz = fz(1,1,:);
rrz = fz(2,1,:);

p=5;
q=5;

% Plot IRF of closed-form solution
PERIODS = 100;
c2 = varma.fit(z,ccz,p,q);
r2 = varma.fit(z,rrz,p,q);
imp = zeros(1,PERIODS);
imp(1) = 1;
res_c = irf(c2,imp,1e-4);
res_r = irf(r2,imp,1e-4);

figure 
subplot(1,2,1)
c0_imp = plot(0:PERIODS-1,res_c,'-b','linewidth',1.2);
 title('IRF of $c_{t}$','Interpreter','latex','Fontsize',12);
set(c0_imp,'LineWidth',2);
xlabel('Periods','Interpreter','latex','Fontsize',12);
ylabel('$(\%)$','Interpreter','latex','Rotation',0);
legend('PI Equilibrium 2, Closed-form Solution','Interpreter','latex');

subplot(1,2,2)
r0_imp = plot(0:PERIODS-1,res_r,'-r','linewidth',1.2);
title('IRF of $r_{kt}$','Interpreter','latex','Fontsize',12);
set(r0_imp,'LineWidth',2);
xlabel('Periods','Interpreter','latex','Fontsize',12);
ylabel('$(\%)$','Interpreter','latex','Rotation',0);
legend('PI Equilibrium 2, Closed-form Solution','Interpreter','latex');
