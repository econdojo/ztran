swaq% This code solves the model of Graham and Wright (2010, JME)  
% Incomplete Information Version Using two equation system: 
% c_{t} -E_{it} [c_{t+1} - (1-beta*(1-delta))*r_{kt+1}] = 0
% r_{kt} = (Pa(z)/Q(z))*a_t + (Pc(z)/Q(z))*c_t
% Solve for Equilibrium 2 where individual consumption is stationary, no
% random walk, idiosyncratic capital has explosive roots at z= beta
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

%% Incomplete Information Iteration
nx = 1;
ne = 2;

% set up intial value
z = linspace(-0.99,0.99,50); % Choosing 100 grid points works as well
nz = length(z);

fz = zeros(1,1,nz);  % restricted analytic function f(z)
Ciz = zeros(nx,1,nz); % idiosyncratic consumption

%% Alternative Initial Value Choices
% for k = 1:nz
%     fz(:,:,k) = 1/(1-rho_a*z(k));
% end  
% for k = 1:nz
%     Ciz(:,:,k) = 1;%1/(1-rho_a*z(k));
% end

%% Iteration Setup

AllIts = 0; % be careful about the caption of iteration count
AllDiff = 1;
ITERATE_ACCURACY = 1e-6;

% exogenous parameters 
r = roots(Qz); % Non-invertible root in the non-expectational block
M = (polyval(Paz,r)*(1/(1-rho_a*r)))/(polyval(Pcz,r)); % Aggregate consumption C_a(r) =M, paper notation S
theta = (1-beta*(1-delta)); % composite parameter in Euler Equation
V = [sig_a^2,0;0,sig_z^2];

% VARMA fitting order
p = 5;
q = 5;

ITERATION_TICK = 10; % show convergence every 100 iterations

tic
while AllDiff >ITERATE_ACCURACY
    
    last_fz = fz;
    last_Ciz = Ciz;
    
    Cz = zeros(nx,1,nz);  % aggregate consumption in terms of aggregate innovations only
    Rkz = zeros(1,ne,nz);  % aggregate capital return in terms of both aggregate and idiosyncratic innovations
    
    for k = 1:nz
        Cz(:,:,k) = 1*(z(k)-r)*fz(:,:,k)+M;  % C_{a} (z) = a (z-r) f(z) +M, restricted subset
    end
    
    for k = 1:nz
        Rkz(:,1,k) = (polyval(Paz,z(k))*(1/(1-rho_a*z(k)))-polyval(Pcz,z(k))*Cz(:,:,k))/(polyval(Qz,z(k))); % aggregate component of rental rate
    end
    
   % Fit VARMA representation for the signal under incomplete information
   Sig = zeros(2,ne,nz);
   
   for k = 1:nz
       Sig(1,1,k) = 1/(1-rho_a*z(k));
       Sig(1,2,k) = 1/(1-rho_z*z(k));
       Sig(2,1,k) = Rkz(1,1,k);
   end
   
   Sig_arma = varma.fit(z,Sig,p,q,true);
   
    % Compute expectation
    predox_Cz = zeros(nx,ne,nz);
    predox_Cz(1,1,:) = Cz; % aggregate component
    predox_Cz(1,2,:) = Ciz;
    
    preCz_arma = varma.fit(z,predox_Cz,p,q,true);
    [Ecz,~] = wh(z,1,Sig_arma, preCz_arma,V,350);  % Expand DFT around 550 points to ensure accuracy
        
    Rkz_arma = varma.fit(z,Rkz,p,q,true);
    [Erz,~] = wh(z,1,Sig_arma,Rkz_arma,V,350);
    
    EEz = zeros(nx,ne,nz);  % Euler equation error
    
    for k = 1:nz
        EEz(:,:,k) = Cz(:,:,k)-Ecz(:,:,k)+theta*Erz(:,:,k);
    end
    
    % Updating 
    fz_temp = zeros(nx,1,nz);  % temporary f(z) function
    Ciz_temp = zeros(nx,1,nz); % temporary idiosyncratic consumption
    
    for k = 1:nz
        fz_temp(:,:,k) = last_fz(:,:,k)+EEz(:,1,k); % aggregate  EE error
        Ciz_temp(:,:,k) = Ecz(1,2,k)-theta*Erz(1,2,k);
        
    end
   
    % set up differential stepsize for each function to improve stability
    l_f = 0.5;
    l_c = 1;
    fz = fz_temp*l_f+last_fz*(1-l_f);
    Ciz = Ciz_temp*l_c+last_Ciz*(1-l_c);
    
     % Compute convergence criteria
     N_1 = max(max(max(abs(last_Ciz))));
     N_2 = max(max(max(abs(last_fz))));
     Norm =[N_1, N_2];
     Diff_1 = max(max(max(abs(fz-last_fz))));
     Diff_2 = max(max(max(abs(Ciz-last_Ciz))));
     AllDiff = max([Diff_1,Diff_2]);
     AllDiff = AllDiff/(max(Norm)+1e-5);
   
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
c2 = varma.fit(z,Cz,p,q,true);
r2 = varma.fit(z,Rkz(:,1,:),p,q,true);
imp = zeros(1,PERIODS);
imp(1) = 1;
res_c = irf(c2,imp);
res_r = irf(r2,imp);

figure 
subplot(1,2,1)
c0_imp = plot(0:PERIODS-1,res_c,'-b','linewidth',1.2);
 title('IRF of $c_{t}$','Interpreter','latex','Fontsize',12);
set(c0_imp,'LineWidth',2);
xlabel('Periods','Interpreter','latex','Fontsize',12);
ylabel('$(\%)$','Interpreter','latex','Rotation',0);
legend('PI Equilibrium 2, Algorithm 3','Interpreter','latex');

subplot(1,2,2)
r0_imp = plot(0:PERIODS-1,res_r,'-r','linewidth',1.2);
title('IRF of $r_{kt}$','Interpreter','latex','Fontsize',12);
set(r0_imp,'LineWidth',2);
xlabel('Periods','Interpreter','latex','Fontsize',12);
ylabel('$(\%)$','Interpreter','latex','Rotation',0);
legend('PI Equilibrium 2, Algorithm 3','Interpreter','latex');
