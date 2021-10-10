% This code solves the model of Graham and Wright (2010, JME)  
% Incomplete Information Version Using two equation system: 
% c_{t} -E_{it} [c_{t+1} - (1-beta*(1-delta))*r_{kt+1}] = 0
% r_{kt} = (Pa(z)/Q(z))*a_t + (Pc(z)/Q(z))*c_t
% Solve for Equilibrium 1 where idiosyncratic consumption and capital features random walk
% Use Algorithm 1 based on first-differenced signal representation
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
z = [-0.999, linspace(-0.95,0.95,100), 0.999] ; % Have to avoid z =beta 

%% Alternative State-Space Choice
%z = [linspace(-0.99,0.99,30)];
%z = [linspace(-0.99,0.93, 68), linspace(0.96,0.999,2)];
%z=linspace(-0.95,0.95,50);

%% Initial Values
nz = length(z); 
Rkz = zeros(nx,ne,nz);

for k = 1:nz
    Rkz(:,1,k) =  0; 
    % Alternative Choices of Initial Values
    %1;    1/(1-0.9*z(k));     %1/(1-rho_a*z(k)); %-1/(1-rho_a*z(k)); %-(z(k)-rho_a);  
end

%% Set up Iterations
AllIts = 0; % be careful about the caption of iteration count
AllDiff = 1;
ITERATE_ACCURACY = 1e-6;

% exogenous parameters 
theta = (1-beta*(1-delta)); % composite parameter in Euler Equation

fa_beta = ((1-beta)*polyval(Paz,beta)/polyval(Pcz,beta)); % Predetermined Constant values of f_az(beta)
fb_beta = -((1-beta)*polyval(Qz,beta))/polyval(Pcz,beta); % Predetermined Constant values of f_bz(beta)

V = [sig_a^2,0;0,sig_z^2];

% VARMA fitting order
p = 5;
q = 5;

ITERATION_TICK = 10; % show convergence every 100 iterations

tic
while AllDiff >ITERATE_ACCURACY
    
   last_Rkz = Rkz;
   
   % Fit VARMA representation for the signal under incomplete information
   Sig = zeros(2,ne,nz);
   
   for k = 1:nz
       Sig(1,1,k) = 1/(1-rho_a*z(k));
       Sig(1,2,k) = 1/(1-rho_z*z(k));
       Sig(2,1,k) = Rkz(1,1,k);
   end
   
   Sig_arma = varma.fit(z,Sig,p,q,true);
   
   % Find Wold representation using state space method, evaluated at z =
   % beta
    SSR = ss(Sig_arma,V,[1 2]);
    ns = size(SSR.F,1);
    
    Gammaz = zeros(2,2,nz);
    
    for k = 1:nz
        Gammaz(:,:,k) = pinv(eye(2)+SSR.H*pinv(eye(ns)-SSR.F*z(k))*SSR.K*z(k));
    end
    
    Gammaz_beta = zeros(2,2,1);
    Gammaz_beta(:,:,1) = pinv(eye(2)+SSR.H*pinv(eye(ns)-SSR.F*beta)*SSR.K*beta);
    
    % Compute Signal Representation of Expectation and Solve for Constants
    Rkz_arma = varma.fit(z,Rkz,p,q,true);
    [~,Erz] = wh(z,1,Sig_arma,Rkz_arma,V,400); %Expand DFT around 400 points to ensure accuracy
    [~,Rke_beta] =  wh(beta,1,Sig_arma,Rkz_arma,V,400);  %Evaluate expectation of r_{kt+1} at z =beta
    
    % First-differenced signal representation used as intermediate step
    faz = zeros(nx,1,nz); % first-differenced analytic functions with respect to the first signal
    fbz = zeros(nx,1,nz); % first-differenced analytic functions with respect to the second signal (endogenous aggregate return of capital)
    
    % Solving for free constants
    G = Gammaz_beta(:,:,1);
    Theta = pinv(G')*[fa_beta-theta*Rke_beta(1,1)*beta;
                      fb_beta-theta*Rke_beta(1,2)*beta];
    
    % Compute first-differenced signal representation function
    for k = 1:nz
        faz(:,:,k) = Theta'*Gammaz(:,1,k)+theta*Erz(:,1,k)*z(k);
        fbz(:,:,k) = Theta'*Gammaz(:,2,k)+theta*Erz(:,2,k)*z(k);
    end

    % Compute Innovation representation for aggregate consumption
      Cz = zeros(nx,1,nz);  % aggregate consumption

      for k = 1:nz
        Cz(:,:,k) = (faz(:,:,k)*polyval(Qz,z(k))+fbz(:,:,k)*polyval(Paz,z(k)))/((1-rho_a*z(k))*(polyval(Qz,z(k))*(1-z(k))+polyval(Pcz,z(k))*fbz(:,:,k)));
      end
      
     % Updating 
      Rkz_temp = zeros(1,ne,nz); % Temporary function for aggregate return to capital 
     
      % Updated Innovation Representation
     for k = 1:nz
        Rkz_temp(:,1,k) = ((polyval(Paz,z(k))*(1-z(k))-polyval(Pcz,z(k))*faz(:,:,k))/((1-rho_a*z(k))*(polyval(Qz,z(k))*(1-z(k))+polyval(Pcz,z(k))*fbz(:,:,k))));
     end
    
    % set up stepsize  to improve stability
    l = 0.3;
    Rkz  = Rkz_temp*l+last_Rkz*(1-l);
    
    % Compute convergence criteria
   Norm = max(max(max(abs(last_Rkz))));
   AllDiff = max(max(max(abs(Rkz-last_Rkz))));
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
legend('PI Equilibrium 1, Algorithm 1','Interpreter','latex');

subplot(1,2,2)
r0_imp = plot(0:PERIODS-1,res_r,'-r','linewidth',1.2);
title('IRF of $r_{kt}$','Interpreter','latex','Fontsize',12);
set(r0_imp,'LineWidth',2);
xlabel('Periods','Interpreter','latex','Fontsize',12);
ylabel('$(\%)$','Interpreter','latex','Rotation',0);
legend('PI Equilibrium 1, Algorithm 1','Interpreter','latex');