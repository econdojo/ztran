% This code solves the Hank model of Angeletos and Huo (2021, AER)  
% Incomplete Information combined with heterogeneity in MPC and Business Cycle Exposure 
% Focus on Two Types of Equilibrium:
% 1. Non-degenerate wealth dynamics, no fiscal transfer
% 2. Network model (no wealth dynamics, with fiscal transfer)
% Allow for Two Information Set
% 1. Exogenous, Symmetric Incomplete Information: Privation signal structure is identical across groups
% 2. Full Information
% In execution, when change parameters please compute the full-info model first, to get
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

% Composite Parameters
beta_1 = beta*w_1;  % Group 1 effective discount factor
beta_2 = beta*w_2;  % Group 2 effective discount factor

theta = pi_2*phi_2*w_1+pi_1*phi_1*w_2; % AR component of fixed point for aggregate consumption (output)


%% Incomplete Information Iteration
nx = 1;
ne = 2;

% set up intial value
z = linspace(-0.99,0.99,100); % Choosing 100 grid points works as well
nz = length(z);

Cz = zeros(1,1,nz);  % initial value for aggregate consumption function

for k = 1:nz
    Cz(:,:,k) =0 ; %(z(k)-0.5)/(1-0.97*z(k));
end

%% Iteration Setup
AllIts = 0; % be careful about the caption of iteration count
AllDiff = 1;
ITERATE_ACCURACY = 1e-6;

V = [sig_r^2,0;0,sig_eta^2];

% VARMA fitting order
p = 5;
q = 5;

% Exogenous Shock 
rz = zeros(nx,ne,nz);

for k = 1:nz
    rz(:,1,k) = 1/(1-rho_r*z(k));
end

rz_arma = varma.fit(z,rz,p,q,true);

% Model and Information set options

wealthdynamics = 0; % If network model (no wealth dynamics with fiscal transfer), set to 0
infor = 1; % symmetric, exogenous incomplete information, full information set to 0


ITERATION_TICK = 10; % show convergence every 100 iterations

tic
while AllDiff >ITERATE_ACCURACY
    
    
    last_Cz = Cz;
    
    
   if infor == 1 
   %Fit VARMA representation for the signal under incomplete information
   Sig = zeros(1,ne,nz);
   
   for k = 1:nz
       Sig(1,1,k) = 1/(1-rho_r*z(k));
       Sig(1,2,k) = 1;
   end
   
   Sig_arma = varma.fit(z,Sig,p,q,true);
   end
   
   
    % Compute expectation
    
    predox_Cz = zeros(nx,ne,nz);
    predox_Cz(1,1,:) = Cz; % aggregate component
    
    preCz_arma = varma.fit(z,predox_Cz,p,q,true);
   

    Ecz_1 = zeros(nx,ne,nz); % Group 1 expectations of aggregate consumption (output), no need to multiply mass here
   
   % Expectation of Discounted future sum with effective discount rate beta_1
   if infor == 1
     [Ecz_1,~] = wh(z,beta_1,Sig_arma, preCz_arma,V,550);  % Expand DFT around 550 points to ensure accuracy
   end
   if infor == 0
   [Ecz_1,~] = wh(z,beta_1,varma({zeros(2)},{eye(2)}), preCz_arma,V,550);
   end
   Ecz_1 = Ecz_1*(1-beta_1)*phi_1;
   
   Ecz_2 = zeros(nx,ne,nz); % Group 1 expectations of aggregate consumption (output), no need to multiply mass here
   
   % Expectation of Discounted future sum with effective discount rate beta_2
   if infor == 1
     [Ecz_2,~] = wh(z,beta_2,Sig_arma, preCz_arma,V,550); % Expand DFT around 550 points to ensure accuracy
   end
   
   if infor == 0
    [Ecz_2,~] = wh(z,beta_2,varma({zeros(2)},{eye(2)}), preCz_arma,V,550);
   end
   Ecz_2 = Ecz_2*(1-beta_2)*phi_2;
           
   Erz_1 = zeros(nx,ne,nz);
   if infor == 1
    [Erz_1,~] = wh(z,beta_1,Sig_arma,rz_arma,V,550);
   end
   
   if infor == 0
     [Erz_1,~] = wh(z,beta_1,varma({zeros(2)},{eye(2)}),rz_arma,V,550);
   end
   
   Erz_1 = Erz_1*beta_1;
   
   Erz_2 = zeros(nx,ne,nz);
   if infor == 1
    [Erz_2,~] = wh(z,beta_2,Sig_arma,rz_arma,V,550);
   end
   
   if infor == 0
   [Erz_2,~] = wh(z,beta_2,varma({zeros(2)},{eye(2)}),rz_arma,V,550);
   end
   Erz_2 = Erz_2*beta_2;
   
   
    % Updating 
    Cz_temp = zeros(1,1,nz); % temporary aggregate consumption
    
    for k = 1:nz
        if wealthdynamics == 1
        Cz_temp(:,:,k) = (1/(1-theta*z(k)))*(pi_1*(1-w_2*z(k))*(Ecz_1(1,1,k)-Erz_1(1,1,k))+pi_2*(1-w_1*z(k))*(Ecz_2(1,1,k)-Erz_2(1,1,k)));
        end
        if wealthdynamics ==0
        Cz_temp(:,:,k) = (pi_1*(Ecz_1(1,1,k)-Erz_1(1,1,k))+pi_2*(Ecz_2(1,1,k)-Erz_2(1,1,k)));
        end
    end
   
    % set up differential stepsize for each function to improve stability
    l = 0.5;
    Cz = Cz_temp*l+last_Cz*(1-l);
    
    % Compute the group wealth dynamics function value
   if wealthdynamics == 1
   S_1z = zeros(nx,1,nz); % group 1 wealth
   S_2z = zeros(nx,1,nz); % group 2 wealth
   
   for k = 1:nz
       S_1z(:,:,k) = (1/(1-w_1*z(k)))*(phi_1*Cz(:,:,k)+Erz_1(1,1,k)-Ecz_1(1,1,k));
       S_2z(:,:,k) = -(pi_1/pi_2)*S_1z(:,:,k);
   end
   end
    
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
PERIODS = 50;
c2 = varma.fit(z,Cz,p,q,true);
imp = zeros(1,PERIODS);
imp(1) = -1;
res_c = irf(c2,imp);

% Full information initial impact for normalization in figure plot, if
% change parameters, always run full-information equilibrium first
if infor == 0 && wealthdynamics ==0
 C_full_impact_1 = res_c(1);
 save('normaliz_1','C_full_impact_1');
end

if infor ==0 && wealthdynamics ==1
 C_full_impact_2 = res_c(1);
  save('normaliz_2','C_full_impact_2');
end


if wealthdynamics == 0
    load('normaliz_1');
    res_c = res_c/C_full_impact_1;
end

if wealthdynamics == 1
    load('normaliz_2');
    res_c = res_c/C_full_impact_2;
end

if wealthdynamics == 1
s1 = varma.fit(z,S_1z,p,q,true);
res_s1 = irf(s1,imp)/C_full_impact_2;
s2 = varma.fit(z,S_2z,p,q,true);
res_s2 = irf(s2,imp)/C_full_impact_2;
end
    

if wealthdynamics == 0
figure 
%subplot(1,2,1)
c0_imp = plot(0:PERIODS-1,res_c,'-b','linewidth',1.2);
 title('IRF of $c_{t}$ ($y_t$)','Interpreter','latex','Fontsize',12);
set(c0_imp,'LineWidth',2);
xlabel('Periods','Interpreter','latex','Fontsize',12);
ylabel('$(\%)$','Interpreter','latex','Rotation',0);
if infor == 1
legend('Hank-Network Model, Symmetric (Exogenous) Incomplete Information','Interpreter','latex');
else
legend('Hank-Network Model, Full Information','Interpreter','latex');
end
end

if wealthdynamics == 1
    
figure 
subplot(1,2,1)
c0_imp = plot(0:PERIODS-1,res_c,'-b','linewidth',1.2);
 title('IRF of $c_{t}$ ($y_t$)','Interpreter','latex','Fontsize',12);
set(c0_imp,'LineWidth',2);
xlabel('Periods','Interpreter','latex','Fontsize',12);
ylabel('$(\%)$','Interpreter','latex','Rotation',0);
if infor == 1
legend('PI Hank+Endogenous Wealth Dynamics','Interpreter','latex');
else
legend('FI Hank+Endogenous Wealth Dynamics','Interpreter','latex');
end


subplot(1,2,2)
r0_imp = plot(0:PERIODS-1,res_s1,'-r',0:PERIODS-1,res_s2,'--g','linewidth',1.2);
title('Endogenous Wealth Dynamics','Interpreter','latex','Fontsize',12);
set(r0_imp,'LineWidth',2);
xlabel('Periods','Interpreter','latex','Fontsize',12);
ylabel('$(\%)$','Interpreter','latex','Rotation',0);
if infor == 1
legend('PI Group 1 Wealth','PI Group 2 Wealth','Interpreter','latex');
else
legend('FI Group 1 Wealth','FI Group 2 Wealth','Interpreter','latex');
end
end

