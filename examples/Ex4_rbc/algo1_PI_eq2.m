% This code solves the model of Graham and Wright (2010, JME)  
% Incomplete Information Version Using two equation system: 
% c_{t} -E_{it} [c_{t+1} - (1-beta*(1-delta))*r_{kt+1}] = 0
% r_{kt} = (Pa(z)/Q(z))*a_t + (Pc(z)/Q(z))*c_t
% Solve for Equilibrium 2 where individual consumption is stationary, no
% random walk, idiosyncratic capital has explosive roots at z= beta
% Use Algorithm 1 based on signal representation
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
z = linspace(-0.99,0.99,100);

%% Alternative State-Space Choice
%z = [linspace(-0.99,0.93, 98), linspace(0.95,0.999,2)];
%z = [-0.99, linspace(-0.95,0.95,98), 0.995] ;
%z =  [linspace(-0.99,0.97,100), 0.999];

%% Initial Values
nz = length(z); % dimension of state space grid
  
faz = zeros(1,1,nz);  % analytic functions with respect to the first signal
for k = 1:nz
    faz(:,:,k) = 0;
    %-(1)/(1-rho_a*z(k)); %1/(1-rho_a*z(k)); %((1-1.2*z(k)))/(1-rho_a*z(k));  %1/(1-rho_a*z(k)); Alternative Choices of Initial values 
end
 
fbz = zeros(nx,1,nz); % analytic functions with respect to the second signal (endogenous aggregate return of capital
  
for k = 1:nz
    fbz(:,:,k) = -z(k); 
    % Alternative Choices of Initial Values
    %-z(k);      % -z(k); %-1;   %-(z(k)+1.2); %-1; %(z(k)-2);
    %-1/(1-rho_a*z(k)); %-z(k); %((1-1.5*z(k)))/(1-rho_a*z(k));  z(k)-0.2;

end

%% Set up Iterations
AllIts = 0; % be careful about the caption of iteration count
AllDiff = 1;
ITERATE_ACCURACY = 1e-6;

theta = (1-beta*(1-delta)); % composite parameter in Euler Equation
V = [sig_a^2,0;0,sig_z^2];

% VARMA fitting order
p = 5; 
q = 5;

ITERATION_TICK = 10; % show convergence every 100 iterations

tic
while AllDiff >ITERATE_ACCURACY
      last_faz = faz;
      last_fbz = fbz;
    
    Cz = zeros(nx,1,nz);  % innovation analytic function of aggregate consumption
    Ciz = zeros(nx,1,nz); % innovation analytic function of idiosyncratic consumption
    Rkz = zeros(1,ne,nz); % innovation analytic function of aggregate capital return in terms of both aggregate and idiosyncratic innovations
    
    % Solve for Innovation Representation Based on Signal Function
    for k = 1:nz
        Cz(:,:,k) = (faz(:,:,k)*polyval(Qz,z(k))+fbz(:,:,k)*polyval(Paz,z(k)))/((1-rho_a*z(k))*(polyval(Qz,z(k))+polyval(Pcz,z(k))*fbz(:,:,k)));
        Ciz(:,:,k) =faz(:,:,k)/(1-rho_z*z(k));
    end
    
    for k = 1:nz
        Rkz(:,1,k) = (polyval(Paz,z(k))-polyval(Pcz,z(k))*faz(:,:,k))/((1-rho_a*z(k))*(polyval(Qz,z(k))+polyval(Pcz,z(k))*fbz(:,:,k))); % Aggregate component only
    end

   % Fit VARMA representation for the signal under incomplete information
   Sig = zeros(2,ne,nz);
   
   for k = 1:nz
       Sig(1,1,k) = 1/(1-rho_a*z(k));
       Sig(1,2,k) = 1/(1-rho_z*z(k));
       Sig(2,1,k) = Rkz(1,1,k);
   end
   
   Sig_arma = varma.fit(z,Sig,p,q);
   
   % Compute signal representation of expectation
    
    predox_Cz = zeros(nx,ne,nz);
    predox_Cz(1,1,:) = Cz; % aggregate component
    predox_Cz(1,2,:) = Ciz;
    
    preCz_arma = varma.fit(z,predox_Cz,p,q);
    [~,Ecz] = wh(z,1,Sig_arma, preCz_arma,V,550); % Expand DFT around 550 points to ensure accuracy

    Rkz_arma = varma.fit(z,Rkz,p,q);
    [~,Erz] = wh(z,1,Sig_arma,Rkz_arma,V,550); 
    
    % Updating 
    faz_temp = zeros(nx,1,nz); % Temporary updating function for faz
    fbz_temp = zeros(nx,1,nz); % Temporary updating function for fbz
    
    for k = 1:nz
        faz_temp(:,:,k) = Ecz(:,1,k)-theta*Erz(:,1,k);
        fbz_temp(:,:,k) = Ecz(1,2,k)-theta*Erz(1,2,k);
    end

    % set up differential stepsize for each function to improve stability
    l_a = 1;
    l_b = 0.7;
    faz = faz_temp*l_a+last_faz*(1-l_a);
    fbz = fbz_temp*l_b+last_fbz*(1-l_b);
    
    % Compute convergence criteria
    N_1 = max(max(max(abs(last_faz))));
    N_2 = max(max(max(abs(last_fbz))));
    Norm =[N_1, N_2];
    Diff_1 = max(max(max(abs(faz-last_faz))));
    Diff_2 = max(max(max(abs(fbz-last_fbz))));
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
legend('PI Equilibrium 2, Algorithm 1','Interpreter','latex');

subplot(1,2,2)
r0_imp = plot(0:PERIODS-1,res_r,'-r','linewidth',1.2);
title('IRF of $r_{kt}$','Interpreter','latex','Fontsize',12);
set(r0_imp,'LineWidth',2);
xlabel('Periods','Interpreter','latex','Fontsize',12);
ylabel('$(\%)$','Interpreter','latex','Rotation',0);
legend('PI Equilibrium 2, Algorithm 1','Interpreter','latex');
