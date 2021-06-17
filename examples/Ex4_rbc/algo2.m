% This code solves the model of Graham and Wright (2010, JME)  
% Using two equation system: 
% c_{t} -E_{it} [c_{t+1} - (1-beta*(1-delta))*r_{kt+1}] = 0
% r_{kt} = (Pa(z)/Q(z))*a_t + (Pc(z)/Q(z))*c_t
% Use Algorithm 2 based on non-expectation iteration with linear
% free-constatnt system (associated with forecast errors)
% the code handles three equilibrium: FI equilibrium, PI equilibrium 1
% (random walk), PI equilibrium 2
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
z = [linspace(-0.99,0.93, 98), linspace(0.97,0.99,2)]; % avoid choices closed to z = r_2

%% Alternative State-Space Choice
%z = linspace(-0.99, 0.999,70);
%z = [-0.999, linspace(-0.92,0.92,100), 0.999] ;

%% Initial Values
nz = length(z);
Rkz = zeros(nx,ne,nz);

for k = 1:nz
    Rkz(:,1,k) = z(k)-r_2; 
    
    % Extensive test of other initial value choices
    %-(z(k)+0.1); %-(1-0.2*z(k)); %-(z(k)-0.1);   
    %-1.3*z(k)+1;  %(1-1.5*z(k));  %z(k)-5; %-1/(1-0.2*z(k)); %(-1-0.8*z(k))*(z(k)+0.9);       
    % (z(k)-1.2)*(z(k)-0.9)/(1-rho_a*z(k)); %-(z(k)-0.2)/(1-rho_a*z(k));     %-(z(k)-0.5);     
    %-1;   %-1/(1-rho_a*z(k)); % (z(k)-1.2)*(z(k)-0.7)/(1-0.9*z(k));  %-(z(k)+0.1);    %(z(k)-r_2)/(1-rho_a*z(k));  
    % -(z(k)+0.9); %-1/(1-rho_a*z(k)) ;       
    %-(z(k)-0.7)/(1-rho_a*z(k)); %/(1-rho_a*z(k)); %        (z(k)-r_2)/(1-rho_a*z(k));   
    % (z(k)-rho_a)*(z(k)-r_2);    % (z(k)-0.6)*(z(k)-1.2)/(1-0.9*z(k));    % -(z(k)-0.5)/(1-0.9*z(k));    %(z(k)-1.5);  %-1/(1-0.8*z(k));
    %(z(k)-0.6)*(z(k)-0.92)/(1-rho_a*z(k))/(1-0.9*z(k));    %(z(k)-0.4)*(z(k)-1.02); %/(1-rho_a*z(k)); 
    %(z(k)-1.2)/(1-rho_a*z(k)) ; (z(k)-1.2) %(z(k)-0.4)    
    %(Rf*(z(i)-r_4)*(z(k)-r_2))/(Phi_1*(1-rho_a*z(i))*(z(i)-r_1));      
    %(z(k)-r_2)*(z(k)-0.9);         %0.4*(z(k)-0.9337)*(z(k)-1.02); 
    %/(1-rho_a*z(k))/(1-0.92*z(k)); %*(z(k)-1.02)
end

%% Set up Iterations

AllIts = 0; % be careful about the caption of iteration count
AllDiff = 1;
ITERATE_ACCURACY = 1e-6;

% exogenous parameters 
theta = (1-beta*(1-delta)); % composite parameter in Euler Equation

V = [sig_a^2,0;0,sig_z^2];

% Predetermined RHS value of linear free constant system associated with
% forecast error
G_beta = -(polyval(Paz,beta)*(1-beta))/(polyval(Pcz,beta)*(1-rho_z*beta));
G_1 = -(polyval(Paz,1)*(1-1))/(polyval(Pcz,1)*(1-rho_z*1));
G_r2 = (theta*polyval(Paz,r_2))/((1-rho_a*r_2)*polyval(Qz,r_2));

% VARMA fitting order
p = 5;
q = 5;

eq_id = 1; % eq_id = 0 (FI equilibrium), 1(PI equilibrium 1), 2 (PI equilibrium 2)

ITERATION_TICK = 100; % show convergence every 100 iterations

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
   
   Sig_arma = varma.fit(z,Sig,p,q);
   
  % Find Wold representation using state space method, and compute
  % endogenous function F(z)
    SSR = ss(Sig_arma,V,[1 2]);
    ns = size(SSR.F,1);

    Fz = zeros(2,2,nz);
    Gammaz = zeros(2,2,nz);
  
  for k = 1:nz
      if eq_id == 0
      Fz(:,:,k) = eye(2);
      else
      Fz(:,:,k) = pinv(eval(0,Sig_arma))*pinv(eye(2)+SSR.H*pinv(eye(ns)-SSR.F*z(k))*SSR.K*z(k))*Sig(:,:,k); %pinv(armaeval(0,Sig_ar,Sig_ma))*
      end
      Gammaz(:,:,k) = eye(2)+SSR.H/(eye(ns)-SSR.F*z(k))*SSR.K*z(k);
  end
  
  Fz_arma = varma.fit(z,Fz,p,q);
  Fz_r2 = eval(r_2,Fz_arma);  % F(z) evaluated at z = r_2
  Fz_beta = eval(beta,Fz_arma); % F(z) evaluated at z = beta
  Fz_1 = eval(1,Fz_arma);  % F(z) evaluated at z = 1
  
  % Solve for free constants
  if eq_id == 1 
  
  F_d = [Fz_r2(1,1), -Fz_r2(2,1);
         Fz_beta(1,2), -Fz_beta(2,2)];
     
  Lambda = pinv(F_d)*[G_r2;
                      G_beta];
  end
  
  if eq_id == 2
      
  F_d = [Fz_r2(1,1), -Fz_r2(2,1);
         Fz_1(1,2), -Fz_1(2,2)];
  Lambda = pinv(F_d)*[G_r2;
                      G_1];
  end
  
  if eq_id == 0
  F_d = eye(2);
  Lambda = pinv(F_d)*[G_r2;
                     G_beta];
  end
 
  % Compute innovation representation for aggregate capital return    
    Rkz_temp = zeros(nx,ne,nz);

    for k = 1:nz
        Rkz_temp(:,1,k) = ((z(k)-1)*polyval(Paz,z(k))/(1-rho_a*z(k))-polyval(Pcz,z(k))*(Lambda(1,1)*Fz(1,1,k)-Lambda(2,1)*Fz(2,1,k)))/polyval(Phiz,z(k)); 
    end

    % set up stepsize  to improve stability
    l = 0.5;
    Rkz = Rkz_temp*l+last_Rkz*(1-l);
     
    % Compute convergence criteria
   N_1 = max(max(max(abs(last_Rkz)))); 
   Norm = N_1;
   AllDiff =  max(max(max(abs(Rkz-last_Rkz))));
   AllDiff = AllDiff/(max(Norm)+1e-5);
   AllIts = AllIts + 1; 
         
         if It_Tick(AllIts)
            Text = ['Inside Iteration: ' num2str(AllIts) ', difference: ' num2str(AllDiff)];
            disp(Text); 
         end   
end
toc
     
% Compute ztran function value for aggregate and idiosyncratic consumption
     Cz = zeros(nx,1,nz);  % aggregate consumption
     Ciz = zeros(nx,1,nz);
     
for k = 1:nz
    Cz(:,1,k) = (polyval(Paz,z(k))-polyval(Qz,z(k))*Rkz_temp(:,1,k)*(1-rho_a*z(k)))/((1-rho_a*z(k))*polyval(Pcz,z(k)));
    Ciz(:,1,k) = (Lambda(1,1)*Fz(1,2,k)-Lambda(2,1)*Fz(2,2,k)); %% first-differenced idiosyncratic consumption to avoid random walk
end
    
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
if eq_id == 1
legend('PI Equilibrium 1, Algorithm 2','Interpreter','latex');
end

if eq_id == 2
legend('PI Equilibrium 2, Algorithm 2','Interpreter','latex');
end

if eq_id == 0
legend('FI Equilibrium, Algorithm 2','Interpreter','latex');
end

subplot(1,2,2)
r0_imp = plot(0:PERIODS-1,res_r,'-r','linewidth',1.2);
title('IRF of $r_{kt}$','Interpreter','latex','Fontsize',12);
set(r0_imp,'LineWidth',2);
xlabel('Periods','Interpreter','latex','Fontsize',12);
ylabel('$(\%)$','Interpreter','latex','Rotation',0);
if eq_id == 1
legend('PI Equilibrium 1, Algorithm 2','Interpreter','latex');
end

if eq_id == 2
legend('PI Equilibrium 2, Algorithm 2','Interpreter','latex');
end

if eq_id == 0
legend('FI Equilibrium, Algorithm 2','Interpreter','latex');
end
