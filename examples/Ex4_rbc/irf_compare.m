%% This code plot the IRF of aggregate consumption in GW(2010) model 
%% Using the frequency domain APFI algorithm, frequency domain closed-form solution
%% and GW(2010) time-domain truncation algorithm based on Nimark (2018)'s method
%% Consider two types of PI equilibrium, and compare accuracy
clear
close all
clc

% Equilibrium 1: Random Walk
% Load IRF sequences from three different methods

load('C_agg_close_eq1');
load('C_agg'); % GW (2010) time-domain sequence
load('C_agg_APFI_eq1');


% Equilibrium 2
% Load IRF sequences from two different methods
load('C_agg_close_eq2');
load('C_agg_APFI_eq2');

PERIODS = 100;
figure 
c0_imp = plot(0:PERIODS-1,C_agg_APFI_eq_1(1,1:PERIODS),'-b',0:PERIODS-1, C_agg_close_eq_1(1,1:PERIODS), 'xm', 0:PERIODS-1, C_agg(1,2:PERIODS+1),'-.g', 'linewidth',1.2);
title('IRF of $c_{t}$: Equilibrium 1','Interpreter','latex','Fontsize',12);
set(c0_imp,'LineWidth',2);
xlabel('Periods','Interpreter','latex','Fontsize',12);
ylabel('$(\%)$','Interpreter','latex','Rotation',0);
legend('APFI algorithm: HTW(2020)','Closed-form solution', 'Time-domain truncation: GW(2010)','Interpreter','latex');


figure
c1_imp = plot(0:PERIODS-1,C_agg_APFI_eq_2(1,1:PERIODS),'-b',0:PERIODS-1, C_agg_close_eq_2(1,1:PERIODS),'xm',  'linewidth',1.2);
title('IRF of $c_{t}$: Equilibrium 1','Interpreter','latex','Fontsize',12);
set(c1_imp,'LineWidth',2);
xlabel('Periods','Interpreter','latex','Fontsize',12);
ylabel('$(\%)$','Interpreter','latex','Rotation',0);
legend('APFI algorithm: HTW(2020)','Closed-form solution','Interpreter','latex');




