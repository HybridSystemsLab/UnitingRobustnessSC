%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: HybridNessterovHBF.m
%--------------------------------------------------------------------------
% Project: Uniting Nesterov's accelerated gradient descent globally with
% heavy ball locally. Nonstrongly convex version.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 09/29/2020 1:07:00

clear all

set(0,'defaultTextInterpreter','tex');

% global variables
global M gamma lambda d beta c_0 c_10 cTilde_0 cTilde_10 d_0 d_10 alpha sigma randomsInterp randomsIndex rho

%%%%%%%%%%%%%%%%%%%%%
%% setting parameters
%%%%%%%%%%%%%%%%%%%%%
setMinima();

% Nesterov constants
kappa = 1; % kappa >= 1, kappa = M/mu = 2/2
M = 2;
d = 1/(sqrt(kappa) + 1);
beta = (sqrt(kappa) - 1)/(sqrt(kappa) + 1);
a = d + beta/(2*kappa);

% Added gradient noise
rho = 0;

% Heavy Ball constants
gamma = 2/3; % 
lambda = 40; 

% Uniting parameters for \mathcal{U}_0 and \mathcal{T}_{1,0}:
c_0 = 3000;  
c_10 = 402.83351565;  

alpha = 1;

% eps_0 has to be bigger than eps_10
eps_0 = 20;
eps_10 = 15;

cTilde_0 = eps_0*alpha
cTilde_10 = eps_10*alpha
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)
d_10 = c_10 - (a^2)*(cTilde_10/alpha)^2 - (1/M)*((cTilde_10^2)/alpha)

%%%%%%%%%%%%%%%%%%%%%
%% Initializing
%%%%%%%%%%%%%%%%%%%%%

% initial conditions
z1_0 = 50; 
z2_0 = 0;
q_0 = 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0];

% simulation horizon
TSPAN=[0 2000]; % 500
TSPAN1=[0 20];
JSPAN = [0 10];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.01);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate the nominal system with same c_10 as noisy solutions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running simulation of nominal system...')

[tNom,jNom,xNom] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

lNom = CalculateLVec(xNom,tNom);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate the perturbed system: sigma = 1 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running noisy simulation with sigma = 1...')

nomLen = 20*length(xNom); % 100*length(xNom);
sample = 0.01; 
   
randoms = randn(1,nomLen);
randInd = 1:nomLen;
randomsInterpIndex = 1:sample:nomLen;
randomsInterp = interp1(randInd,randoms,randomsInterpIndex);
randomsIndex = 1;
sigma = 1; 

[ta,ja,xa] = HyEQsolver(@fUNoise,@gU,@CUNoise,@DUNoise,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

lA = CalculateLVec(xa,ta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simulate the perturbed system: sigma = 5 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running noisy simulation with sigma = 5...')

randomsIndex = 1;
sigma = 5; 

[tb,jb,xb] = HyEQsolver(@fUNoise,@gU,@CUNoise,@DUNoise,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

lB = CalculateLVec(xb,tb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate the perturbed system: sigma = 10 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running noisy simulation with sigma = 10...')

randomsIndex = 1;
sigma = 10; 

[tc,jc,xc] = HyEQsolver(@fUNoise,@gU,@CUNoise,@DUNoise,...
    x0,TSPAN,JSPAN,rule,options,'ode45');
 
lC = CalculateLVec(xc,tc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate the perturbed system: sigma = 15 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running noisy simulation with sigma = 15...')

randomsIndex = 1;
sigma = 15; 

[td,jd,xd] = HyEQsolver(@fUNoise,@gU,@CUNoise,@DUNoise,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

lD = CalculateLVec(xd,td);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate the perturbed system: sigma = 20 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running noisy simulation with sigma = 20...')

randomsIndex = 1;
sigma = 20; 
 
[te,je,xe] = HyEQsolver(@fUNoise,@gU,@CUNoise,@DUNoise,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

lE = CalculateLVec(xe,te);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate the perturbed system: sigma = 25 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running noisy simulation with sigma = 25...')

randomsIndex = 1;
sigma = 25; 

[tf,jf,xf] = HyEQsolver(@fUNoise,@gU,@CUNoise,@DUNoise,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

lF = CalculateLVec(xf,tf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate the perturbed system: sigma = 0.01 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running noisy simulation with sigma = 0.01...')

randomsIndex = 1;
sigma = 0.01; 

[tg,jg,xg] = HyEQsolver(@fUNoise,@gU,@CUNoise,@DUNoise,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

lG = CalculateLVec(xg,tg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate the perturbed system: sigma = 0.1 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running noisy simulation with sigma = 0.1...')
randomsIndex = 1;
sigma = 0.1; 

[th,jh,xh] = HyEQsolver(@fUNoise,@gU,@CUNoise,@DUNoise,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

lH = CalculateLVec(xh,th);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate the perturbed system: sigma = 0.5 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running noisy simulation with sigma = 0.5...')

randomsIndex = 1;
sigma = 0.5; 

[ti,ji,xi] = HyEQsolver(@fUNoise,@gU,@CUNoise,@DUNoise,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

lI = CalculateLVec(xi,ti);