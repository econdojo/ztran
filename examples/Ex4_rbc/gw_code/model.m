clear;
close all;

global r h sigma_c ITERATION_TICK MAX_ITS

% set up system in form
% W(t+1) = G_w W(t) + G_c c(t) + G_y y(t) + eta(t+1)
% M_y y(t) = M_w' W(t) + M_c c(t) + M_om om(t)
% om(t) = Rho_m * om (t-1) + eps(t)

%*******************************************************************************************
%******************************  VARIABLES        ******************************************
%*******************************************************************************************
XiCount = 0;
%states - aggregates before idiosyncratic
XiCount = XiCount + 1;
i_k =XiCount;
XiNames='Aggregate Capital';
XiCount = XiCount + 1;
i_a =XiCount;
XiNames=strvcat(XiNames,'Aggregate Technology');
XiCount = XiCount + 1;
i_ksd =XiCount;
XiNames=strvcat(XiNames,'Household Capital deviation');
XiCount = XiCount + 1;
i_zs =XiCount;
XiNames=strvcat(XiNames,'Idiosyncratic Technology');

% static controls
YCount = 0;
YCount = YCount + 1;
i_r =YCount;
YNames='interest rate';
YCount = YCount + 1;
i_y =YCount;
YNames=strvcat(YNames,'Aggregate Output');
YCount = YCount + 1;
i_w =YCount;
YNames=strvcat(YNames,'Aggregate Wage');
YCount = YCount + 1;
i_l =YCount;
YNames=strvcat(YNames,'Aggregate Employment');
YCount = YCount + 1;
i_wsd =YCount;
YNames=strvcat(YNames,'idio wage deviation');
YCount = YCount + 1;
i_lsd =YCount;
YNames=strvcat(YNames,'idio employment deviation');
YCount = YCount + 1;
i_ks =YCount;
YNames=strvcat(YNames,'Idiosyncratic Capital');
YCount = YCount + 1;
i_ws =YCount;
YNames=strvcat(YNames,'Idiosyncratic Wage');

%  dynamic controls
CCount = 0;
CCount = CCount + 1;
i_c =CCount;
CNames='Aggregate Consumption';
CCount = CCount + 1;
i_cs =CCount;
CNames=strvcat(CNames,'Idio cons');

% specify aggregates
Y_agg = [i_r i_y i_w i_l];
Xi_agg = [i_k i_a];

% specify what's observable
OBSERVABLES = [i_r  i_ws  i_ks ];
MEAS_ERROR = [   i_ws  ];
%*******************************************************************************************
%******************************  CALIBRATION        ******************************************
%*******************************************************************************************
R_bar= 0.015;
g=0; % set syeady-state growth to 0
alpha=.667;
delta=.025;
N_bar = 1/3;
sigma_c = 1;
gamma = 5;

% shocks
phi_agg = 0.9;
var_agg = 0.7^2;

phi_idio= 0.9;
var_idio = 4.9^2;

var_meas = 0.00001;        

h = 5;     

PERIODS = 50;

%*******************************************************************************************
%**********************************  OUTPUT       ******************************************
%*******************************************************************************************
PLOT_PERIODS = 50;

PLOT_CERTAIN = 0; 
PLOT_CERTAIN_DEV = 0;
PLOT_FULL_INF = 0;
PLOT_CAPITAL = 0;
PLOT_LABOUR = 0;
PLOT_HIERACHY = 0;  % plots all the orders for the specified variable eg PLOT_HIERACHY = 2 plots technology
PLOT_FIRST_ORDER = 1;
PLOT_C_IDIO = 0;
PLOT_C_AGG = 1;

SIMULATE = 0;
UNCOND_VAR = 0;     % calculate the unconditional variance of the states
HIERARCHY_IMPACT = 0;   % display difference between solution and certainty equivalent values
%*******************************************************************************************
%**********************************  SETTINGS       ******************************************
%*******************************************************************************************

IT_STYLE = 3;               % 1 = separate iterations, 2 = eta separate, 3 = one iteration
ITERATION_TICK = 100;
MAX_ITS = 20000;
ITERATE_ACCURACY = 1e-5;  
DEBUG  =0;
ImpactCount = 0;
RESTRICTED_ETA = 0;     % this attempts to put some structure on eta to improve convergence.

    
%*******************************************************************************************
%******************************  SMART STARTING VALUES**************************************
%*******************************************************************************************
It_start = 4;
It_end = 1;
It_step = -0.0001;
Range = [It_start:It_step:It_end];
%for var_val =  Range
%var_idio = var_val;
SAVE_VARS = 1;
LOAD_VARS = 0;


%*******************************************************************************************
%******************************  MATRICES        ******************************************
%*******************************************************************************************
G_w = zeros (XiCount,XiCount);
G_c = zeros (XiCount,1);
G_s = zeros (XiCount,1);
G_y = zeros (XiCount,YCount);
M_y = eye(YCount);
M_w= zeros (YCount,XiCount);
M_c= zeros (YCount,1);
M_s= zeros (YCount,1);
M_om = eye(YCount);
Q_w =zeros(XiCount,XiCount);

%*******************************************************************************************
%******************************  PARAMETERS       ******************************************
%*******************************************************************************************
sigma_n = 1/gamma;
squig = sigma_n*(1-N_bar)/N_bar;
nu = squig/(1+(1-alpha)*squig);

% sigma_n = 0 means squig = 0, nu = 0
% sigma_n = inf, squig = inf, nu = 1/(1-alpha)

lambda1= (1+R_bar)/(1+g);
lambda2= alpha*(R_bar+delta)/(1-alpha)/(1+g);
lambda3 =alpha*(R_bar+delta)/(1+R_bar);
lambda4= 1-lambda1-lambda2;

% calibrating theta
K_over_Y = (1-alpha)/(R_bar+delta);
C_over_Y = 1- delta*K_over_Y;
theta_util = 1/(C_over_Y*N_bar*(1-N_bar)^(-sigma_n));

%*******************************************************************************************
%******************************  STATES           ******************************************
%*******************************************************************************************

% aggregate capital - A16
% k(t+1) = lambda1 k(t) + lambda2 (a(t) + n(t)) + lambda4 c(t)
G_w(i_k,i_k)=lambda1;
G_w(i_k,i_a)=lambda2;
G_y(i_k,i_l) = lambda2;
G_c(i_k,1)=lambda4;

% idiosyncratic capital (in deviations) - A17
% ksd(t+1) = lambda1 ksd(t) + lambda2 (zs(t) + nsd(t)) + lambda4 (cs(t)-c(t))
G_w(i_ksd,i_ksd)=lambda1;
G_w(i_ksd,i_zs)=lambda2;
G_y(i_ksd,i_lsd) = lambda2;
G_c(i_ksd,1)=-lambda4;
G_s(i_ksd,1)=lambda4;

% aggregate technology - 12
G_w(i_a,i_a)=phi_agg;
Q_w(i_a,i_a)=var_agg;

% idiosyncratic technology - 13
G_w(i_zs,i_zs)=phi_idio;
Q_w(i_zs,i_zs)=var_idio;

%*******************************************************************************************
%******************************  MEASUREMENT      ******************************************
%*******************************************************************************************
% real interest rate - A9
% r(t) = lambda3 * (a(t) - a(t) + n(t))
M_w(i_r,i_k) = -lambda3;
M_w(i_r,i_a) = lambda3;
M_y(i_r,i_l) = -lambda3; % on LHS

% aggregate labour demand - A10
% w(t) = y(t) - n(t)
M_y(i_w,i_l) = 1; % on LHS
M_y(i_w,i_y) = -1; 

% aggregate output - A11
% y(t) = alpha (a(t) + n(t)) +(1-alpha) k(t)
M_w(i_y,i_k) = 1-alpha;
M_w(i_y,i_a) = alpha;
M_y(i_y,i_l) = -alpha; % on LHS

% this just checks the substitution for the two labour supply equations
% SUBST = 1 uses the substitution, equations A20 and A22
% SUBST = 0 uses the "raw" FOCs, equations A19 and A21

SUBST = 0;

if ~SUBST
    % aggregate labour supply - A19
    %l(t) = squig (w(t) - c(t))
    M_y(i_l,i_w) = -squig;  % on LHS
    M_c(i_l,1) = -squig;

    % household labour supply  - A21      
    %lsd(t) =  nu((wsd(t) - cs(t)+c(t))
    M_y(i_lsd,i_wsd) = -squig; % on LHS
    M_c(i_lsd,1) = squig;
    M_s(i_lsd,1) = -squig;
else
    % aggregate labour supply - A20
    %l(t) = nu ((1-alpha) k(t) + alpha a(t) - c(t))
    M_w(i_l,i_k) = nu * (1-alpha);
    M_w(i_l,i_a) = nu * alpha;
    M_c(i_l,1) = -nu;

    % household labour supply  - A22      
    M_w(i_lsd,i_zs) = squig;
    M_s(i_lsd,1) = -squig;
    M_c(i_lsd,1) = squig;
end

% firm labour demand - A15
% wsd(t) = ysd(t) - lsd(t) = zsd(t)
M_y(i_wsd,i_wsd) = 1; % on LHS
M_w(i_wsd,i_zs) = 1;

% household capital (observable)
% ks = ksd + k
M_y(i_ks,i_ks) = 1; % on LHS
M_w(i_ks,i_k) = 1; 
M_w(i_ks,i_ksd) = 1; 

% household wage (observable)
% ws = wsd + w
M_y(i_ws,i_ws) = 1; % on LHS
M_y(i_ws,i_wsd) = -1; % on LHS
M_y(i_ws,i_w) = -1; % on LHS

%*******************************************************************************************
%******************************  FORECASTING R    ******************************************
% ******************************************************************************************

r_index = i_r;
ETA_CHECK = 0;

Translate
HierachySolver
Analytics

if isnan(eta)
    eta
    error('eta not defined');
end

if isnan(beta)
    beta
    error('beta not defined');
end

Simul

if HIERARCHY_IMPACT
    ImpactCount = ImpactCount + 1;
    Impact = 0;
    for i = 1:PERIODS
        Impact = Impact + (C_agg(i) -  C_certain(i))^2;
    end
    disp('Impact of hierarchy');
    disp(Impact);
    HierarchyImpact(ImpactCount,1) = var_idio;
    HierarchyImpact(ImpactCount,2) = Impact;
end

C_agg(2)

%end % for smart starting values


PlotIRF



if UNCOND_VAR
    %VAR = CalcUncondVar(AR,NN,Q);
    VAR_NEW = CalcUncondVarIter(AR,NN,Q);

    VAR_W = reshape(diag(VAR_NEW),r,2*h+1);
    % exclude idiosync capital
    VAR_W(1:2,1:h+1)
    VAR_W(4,1:h+1)
end



