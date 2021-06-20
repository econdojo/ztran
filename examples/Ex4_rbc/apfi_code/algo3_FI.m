% This code solves the model of Graham and Wright (2010, JME)  
% Using two equation system: 
% c_{t} -E_{it} [c_{t+1} - (1-beta*(1-delta))*r_{kt+1}] = 0
% r_{kt} = (Pa(z)/Q(z))*a_t + (Pc(z)/Q(z))*c_t
% Solve for Full Information Equilibrium (only aggregate part, with
% aggregate TFP shock)
% Use Algorithm 3 based on Euler Equation Error Updating
% No steady-state growth, g =0

clear
close all
clc
global ITERATION_TICK

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

%% Model primitive structure, follow supplementary appendix analysis

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

%% APFI Iteration

nx = 1;
ne = 1;

%set up intial value
z = linspace(-0.99,0.99,50);
nz = length(z);
fz = zeros(nx,ne,nz); %restricted analytic function f(z) associated with the aggregate consumption
for k = 1:nz
    fz(:,:,k) = 0; %1/(1-rho_a*z(k));
end

%% Iteration Setup
AllIts = 0; % be careful about the caption of iteration count
AllDiff = 1;
ITERATE_ACCURACY = 1e-6;

% exogenous parameters 
r = roots(Qz); % Non-invertible root in the non-expectational block
M = (polyval(Paz,r)*(1/(1-rho_a*r)))/(polyval(Pcz,r)); % Aggregate consumption C_a(r) =M, paper notation S
theta = (1-beta*(1-delta)); % composite parameter in Euler Equation

% VARMA fitting order
p = 5;
q = 5;

ITERATION_TICK = 10; % show convergence every 100 iterations

tic
while AllDiff >ITERATE_ACCURACY
    
    last_fz = fz;
    
    Cz = zeros(nx,ne,nz); % Innovation representation of aggregate consumption
    Rkz = zeros(nx,ne,nz); % Innovation representation of aggregate return of capital
    
    for k = 1:nz
        Cz(:,:,k) = 0.5*(z(k)-r)*fz(:,:,k)+M; % C_{a} (z) = a (z-r) f(z) +M, restricted subset, a =0.5
    end
    
    for k = 1:nz
        Rkz(:,:,k) = (polyval(Paz,z(k))*(1/(1-rho_a*z(k)))-polyval(Pcz,z(k))*Cz(:,:,k))/(polyval(Qz,z(k)));
    end
    
    % Compute expectation
    Cz_arma = varma.fit(z,Cz,p,q);
    [Ecz,~] = wh(z,1,varma({0},{1}),Cz_arma,sig_a^2,350);   % Et(ct+1)
    
    Rkz_arma = varma.fit(z,Rkz,p,q);
    [Erz,~] = wh(z,1,varma({0},{1}),Rkz_arma,sig_a^2,350);
    
    EEz = zeros(nx,ne,nz);  % Euler equation error
    
    for k = 1:nz
        EEz(:,:,k) = Cz(:,:,k)-Ecz(:,:,k)+theta*Erz(:,:,k);
    end

    % Updating by EE error
    fz_temp = zeros(nx,ne,nz);
    
    for k = 1:nz
        fz_temp(:,:,k) = last_fz(:,:,k)+EEz(:,:,k);
    end

    % set up  stepsize for each function to improve stability
    l = 0.5;
    fz = fz_temp*l+last_fz*(1-l);
    
    % Compute convergence criteria
    AllDiff = max(max(max(abs(fz-last_fz))));
    AllDiff = AllDiff/(max(max(max(abs(last_fz))))+1e-5);
    AllIts = AllIts + 1; 
         
         if It_Tick(AllIts)
            Text = ['Inside Iteration: ' num2str(AllIts) ', difference: ' num2str(AllDiff)];
            disp(Text); 
         end
    
end
toc

% Plot IRF of Solution

% MA expansion of aggregate consumption and capital return
PERIODS = 100;
c2 = varma.fit(z,Cz,p,q);
r2 = varma.fit(z,Rkz(:,1,:),p,q);
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
legend('FI Equilibrium, Algorithm 3','Interpreter','latex');

subplot(1,2,2)
r0_imp = plot(0:PERIODS-1,res_r,'-r','linewidth',1.2);
title('IRF of $r_{kt}$','Interpreter','latex','Fontsize',12);
set(r0_imp,'LineWidth',2);
xlabel('Periods','Interpreter','latex','Fontsize',12);
ylabel('$(\%)$','Interpreter','latex','Rotation',0);
legend('FI Equilibrium, Algorithm 3','Interpreter','latex');
