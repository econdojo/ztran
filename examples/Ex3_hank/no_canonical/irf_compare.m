%% This code plot the IRF of aggregate consumption and wealth dynamics in AH(2021) AER model 
%% Using the frequency domain APFI algorithm
%% Consider two types of equilibrium: HANK-Network Model with no wealth dynamics, HANK Modl with wealth dynamics
%% Consider two types of information set: Full Information,  Symmetric Incomplete Information with exogenous private signal
clear
close all
clc

% IRF of Aggregate Consumption under Different Model and Information Sets
% Load IRF sequences from four different cases
load('C_FI_nw');
load('C_PI_nw');
load('IRF_FI_w');
load('IRF_PI_w');




PERIODS = 30;
figure 
subplot(1,2,1)
c0_imp = plot(0:PERIODS-1,C_FI_nw(1,1:PERIODS),'-b',0:PERIODS-1, C_FI_w(1,1:PERIODS), '-r', 0:PERIODS-1, C_PI_nw(1,1:PERIODS),'-.g', 0:PERIODS-1, C_PI_w(1,1:PERIODS),'-.k', 'linewidth',1.2);
title('IRF of $c_{t}$ ($y_t$)','Interpreter','latex','Fontsize',12);
set(c0_imp,'LineWidth',2);
xlabel('Periods','Interpreter','latex','Fontsize',12);
ylabel('$(\%)$','Interpreter','latex','Rotation',0);
legend('FI HANK-Network','FI HANK-Endogenous Wealth', 'PI HANK-Network','PI HANK-Endogenous Wealth','Interpreter','latex');

subplot(1,2,2)
s0_imp = plot(0:PERIODS-1,S1_FI_w(1,1:PERIODS),'-b',0:PERIODS-1, S1_PI_w(1,1:PERIODS),'--r', 0:PERIODS-1,S2_FI_w(1,1:PERIODS),'-b', 0:PERIODS-1, S2_PI_w(1,1:PERIODS),'--r', 'linewidth',1.2);
title('Endogenous Wealth Dynamics','Interpreter','latex','Fontsize',12);
set(s0_imp,'LineWidth',2);
xlabel('Periods','Interpreter','latex','Fontsize',12);
ylabel('$(\%)$','Interpreter','latex','Rotation',0);
legend('FI Group 1 Wealth','PI Group 1 Wealth','FI Group 1 Wealth', 'PI Group 2 Wealth','Interpreter','latex');




