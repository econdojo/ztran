% This code solves the Hank model of Angeletos and Huo (2021, AER)  
% Incomplete Information combined with heterogeneity in MPC and Business Cycle Exposure 
% Focus on Equilibrium:  Network model (no wealth dynamics, with fiscal transfer)
% Allow for Two Information Set
% 1. Asymmetric Incomplete Information: Privation signal structure is
% identical across groups, with one group (group 1) recieving a noisy endogenous
% signal
% 2. Full Information
% In execution, when change parameters please compute the full-infor model first, to get
% normalization value for figure plot

clear
close all
clc
global ITERATION_TICK

% Parameters
beta = 0.99;        % quarterly discount rate
m1 = 0.55;          % MPC of Group 1
m2 = 0.05;          % MPC of Grop 2
phi_1 = 2;          % Group 1 Business Cycle exposure
phi_2 = 0;          % Group 2 Business Cycle exposure
pi_1 = 0.5;         % Group 1 mass
pi_2 = 0.5;         % Group 2 mass
w_1 = (1-m1)/beta;  % Group 1 survival probability
w_2 = (1-m2)/beta;  % Group 2 survival probability


% Shock parameters
rho_r = 0.95;        % persistence of aggregate interest rate shock
sig_r = 1;        % std. dev. of aggregate interest rate
sig_eta = 1.61;        % std. dev. of idiosyncratic noise in private signal
sig_m = 4;       % std. dev. of group-specific noise in endogenous signal for group 1
rho_m = 0.95;    % persistence of noise in endogenous signal

% Composite Parameters
beta_1 = beta*w_1;  % Group 1 effective discount factor
beta_2 = beta*w_2;  % Group 2 effective discount factor

theta = pi_2*phi_2*w_1+pi_1*phi_1*w_2; % AR component of fixed point for aggregate consumption (output)


%% Incomplete Information Iteration
nx = 1;
ne = 3;

% set up intial value
z = linspace(-0.99,0.99,50); % Choosing 100 grid points works as well
nz = length(z);

Cz = zeros(nx,ne-1,nz);  % initial value for aggregate consumption function (two shocks)

for k = 1:nz
    Cz(:,:,k) = zeros(nx,ne-1); %(z(k)-0.5)/(1-0.97*z(k));
end


%% Iteration Setup

AllIts = 0; % be careful about the caption of iteration count
AllDiff = 1;
ITERATE_ACCURACY = 1e-6;


V = [sig_r^2,0,0;0,sig_m^2,0;0,0,sig_eta^2];


% VARMA fitting order
p = 5;
q = 5;

% Exogenous Interest Rate Shock 
rz = zeros(nx,ne,nz);

for k = 1:nz
    rz(:,1,k) = 1/(1-rho_r*z(k));
end

rz_arma = varma.fit(z,rz,p,q,true);

% Information set options

%wealthdynamics = 0; % If network model (no wealth dynamics with fiscal transfer), set to 0
infor = 1; % Non-symmetric incomplete information, full information set to 0


ITERATION_TICK = 10; % show convergence every 100 iterations

tic
while AllDiff >ITERATE_ACCURACY
    
    
    last_Cz = Cz;
    
    
   if infor == 1 
   %Fit VARMA representation for the signal under incomplete information
   Sig_1 = zeros(2,ne,nz); % Signal system for group 1
   
   for k = 1:nz
       Sig_1(1,1,k) = 1/(1-rho_r*z(k));
       Sig_1(1,3,k) = 1;
       Sig_1(2,1,k) = Cz(:,1,k);
       Sig_1(2,2,k) = Cz(:,2,k)+1/(1-rho_m*z(k));
   end
   
   Sig1_arma = varma.fit(z,Sig_1,p,q,true);
   
   Sig_2 = zeros(1,ne,nz);
  
   for k = 1:nz
       Sig_2(1,1,k) = 1/(1-rho_r*z(k));
       Sig_2(1,3,k) = 1;
   end
   
   Sig2_arma = varma.fit(z,Sig_2,p,q,true);
   
   end
   
   
    % Compute expectation
    
    predox_Cz = zeros(nx,ne,nz);
    predox_Cz(1,1:2,:) = Cz; % aggregate component
    
    preCz_arma = varma.fit(z,predox_Cz,p,q,true);
   

    Ecz_1 = zeros(nx,ne,nz); % Group 1 expectations of aggregate consumption (output), no need to multiply mass here
   
   % Expectation of Discounted future sum with effective discount rate beta_1
   if infor == 1
     [Ecz_1,~] = wh(z,beta_1,Sig1_arma, preCz_arma,V,350);  
   end
   if infor == 0
   [Ecz_1,~] = wh(z,beta_1,varma({zeros(3)},{eye(3)}), preCz_arma,V,350);
   end
   
   Ecz_1 = Ecz_1*(1-beta_1)*phi_1;
   
   Ecz_2 = zeros(nx,ne,nz); % Group 2 expectations of aggregate consumption (output), no need to multiply mass here
   
   % Expectation of Discounted future sum with effective discount rate beta_2
   if infor == 1
     [Ecz_2,~] = wh(z,beta_2,Sig2_arma, preCz_arma,V,350); % Use 350 points to increase speed
   end
   
   if infor == 0
    [Ecz_2,~] = wh(z,beta_2,varma({zeros(3)},{eye(3)}), preCz_arma,V,350);
   end
   Ecz_2 = Ecz_2*(1-beta_2)*phi_2;
           
   Erz_1 = zeros(nx,ne,nz);
   if infor == 1
    [Erz_1,~] = wh(z,beta_1,Sig1_arma,rz_arma,V,350);
   end
   
   if infor == 0
     [Erz_1,~] = wh(z,beta_1,varma({zeros(3)},{eye(3)}),rz_arma,V,350);
   end
   
   Erz_1 = Erz_1*beta_1;
   
   Erz_2 = zeros(nx,ne,nz);
   if infor == 1
    [Erz_2,~] = wh(z,beta_2,Sig2_arma,rz_arma,V,350);
   end
   
   if infor == 0
   [Erz_2,~] = wh(z,beta_2,varma({zeros(3)},{eye(3)}),rz_arma,V,350);
   end
   Erz_2 = Erz_2*beta_2;
   
   
    % Updating 
    Cz_temp = zeros(nx,ne-1,nz); % temporary aggregate consumption
    
    for k = 1:nz

        Cz_temp(:,:,k) = (pi_1*(Ecz_1(1,1:2,k)-Erz_1(1,1:2,k))+pi_2*(Ecz_2(1,1:2,k)-Erz_2(1,1:2,k)));
        
    end
   
    % set up differential stepsize for each function to improve stability
    l = 1;
    Cz = Cz_temp*l+last_Cz*(1-l);
    

 % Compute convergence criteria
 
 Norm =max(max(max(abs(last_Cz))));
 AllDiff = max(max(max(abs(Cz-last_Cz))));
 AllDiff = AllDiff/(max(Norm)+1e-5);
   
AllIts = AllIts + 1; 
         
         if It_Tick(AllIts)
            Text = ['Inside Iteration: ' num2str(AllIts) ', difference: ' num2str(AllDiff)];
            disp(Text); 
         end   
end
toc


% Plot IRF of Solution
% MA expansion of aggregate consumption 
% Minimum level of IRF on display is set to 1e-4;
PERIODS = 50;
c2 = varma.fit(z,Cz,p,q,true);
imp = zeros(2,PERIODS);
imp(1,1) = -1;
C_agg_1 = irf(c2,imp);
imp = zeros(2,PERIODS);
imp(2,1) = -1;
C_agg_2 = irf(c2,imp);

% Full information initial impact for normalization in figure plot, if
% change parameters, always run full-information equilibrium first
if infor == 0 
 C_full_impact_3 = C_agg_1(1);
 save('normaliz_3','C_full_impact_3');
end

load('normaliz_3');
C_agg_1 = C_agg_1/C_full_impact_3;
C_agg_2 = C_agg_2/C_full_impact_3;

 
if infor == 1
figure 
subplot(1,2,1)
c0_imp = plot(0:PERIODS-1,C_agg_1,'-b','linewidth',1.2);
 title('IRF of $c_{t}$ ($y_t$): Fundamental Shock','Interpreter','latex','Fontsize',12);
set(c0_imp,'LineWidth',2);
xlabel('Periods','Interpreter','latex','Fontsize',12);
ylabel('$(\%)$','Interpreter','latex','Rotation',0);

legend('Non-Symmetric PI Hank-Network Model','Interpreter','latex','Fontsize',7);
%legend('Hank-Network Model, Full Information','Interpreter','latex');

subplot(1,2,2)
c1_imp = plot(0:PERIODS-1,C_agg_2,'-r','linewidth',1.2);
title('IRF of $c_{t}$ ($y_t$): Noise Shock to Group 1','Interpreter','latex','Fontsize',12);
set(c1_imp,'LineWidth',2);
xlabel('Periods','Interpreter','latex','Fontsize',12);
ylabel('$(\%)$','Interpreter','latex','Rotation',0);
% if infor == 1
legend(' Non-Symmetric PI Hank-Network Model','Interpreter','latex','Fontsize',7);
% else
% legend('Hank-Network Model, Full Information','Interpreter','latex');

else 
figure 
c0_imp = plot(0:PERIODS-1,C_agg_1,'-b','linewidth',1.2);
 title('IRF of $c_{t}$ ($y_t$): Fundamental Shock','Interpreter','latex','Fontsize',12);
set(c0_imp,'LineWidth',2);
xlabel('Periods','Interpreter','latex','Fontsize',12);
ylabel('$(\%)$','Interpreter','latex','Rotation',0); 
legend('Hank-Network Model, Full Information','Interpreter','latex','Fontsize',7);
end

