% Optimal discretionary monetary policy in a New Keynesian (NK) model
% with 3 equations using the Taylor rule
% Valerio Dionisi, v.dionisi@campus.unimib.it
% Department of Economics, Management and Statistics - University of Milano-Bicocca

%--------------------------------------------------------------------------

% Situation: the Central Bank (CB) discretionaly adjusts the interest
%            rate to reduce losses using the Taylor rule

%%% 1) PREAMBLE - Parameter Calibration

var x pi i v s;  % endogenous
varexo ss;       % exogenous

% x : output gap (different from y: output)
% pi: rate of inflation
% i : nominal interest rate
% v : monetary policy innovation
% s : markup-shock
% ss: markup-shock innovation

parameters betta sig phin tetap alpha epsi omega kapa rhos phi1 phi2 phi3; % parameters

% betta: discount factor
% sig  : coefficient of relative risk aversion
% tetap: probability of not revise prices in Calvo model
% ...
% phi1 : coefficient of interest rate smoothing
% phi2 : T.R. coefficient on inflation
% phi3 : T.R. coefficient on output gap 

betta = 0.99;
sig   = 1;
phin  = 1;
tetap = 0.75;
alpha = 0.33;
epsi  = 4;
omega = (1-alpha)/(1-alpha+alpha*epsi);
kapa  = ((1-tetap)*(1-tetap*betta)/tetap)*omega*(sig+(phin+alpha)/(1-alpha));

rhos  = 0.6;

% monetary rule
phi1 = 0;       
phi2 = 1.5;
phi3 = 1.5/4;

%% 2) MODEL SETTING  - Write the equations of the model, expressing all at time t
%                      Note the log-linearized model (no steady state and linear specification)

model(linear);

x = x(+1) - 1/sig*(i - pi(+1));  % New-Keynesian IS Curve (NKIS)

pi = betta*pi(+1) + kapa*x + s;  % New-Keynesiam Phillips Curve (NKPC)

s = rhos*s(-1) + ss;             % mark-up shock

% monetary policy instrument
i = phi1*i(-1) + (1-phi1)*(phi2*pi+phi3*x) + v;  % Taylor rule (instead of the Fisher equation)

end;

%% 3) No INITIAL VALUES since I am working with a log-linearized model

%% 4) CHECK STEADY STATE DIFFERENCES - Differences in computed steady state
steady;
resid;

%% 5) SHOCKS 

shocks;
var ss = 0.01^2;
end;

% stoch_simul(irf=40);  % IRFs under no monetary intervention

%% 6) OPTIMAL DISCRETIONARY MONETARY POLICY

planner_objective pi^2 + kapa/epsi*x^2;

discretionary_policy(instruments=(i),irf=40,planner_discount=betta,discretionary_tol=1e-12) x pi i v s;

% This is what happens in the situation in which it is not implemented the Fisher equation, but the Taylor rule

%--------------------------------------------------------------------------
