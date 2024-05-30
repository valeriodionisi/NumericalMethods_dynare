% Optimal commitment monetary policy in a New Keynesian (NK) model
% with 3 equations using the Fisher equation
% Valerio Dionisi, v.dionisi@campus.unimib.it
% Department of Economics, Management and Statistics - University of Milano-Bicocca

%--------------------------------------------------------------------------

% Situation: the Central Bank (CB) is committed(Rmasey planner problem) to adjust the
%            adjust the interest rate to reduce losses using the Fisher equation

%%% 1) PREAMBLE - Parameter Calibration

var x pi i rr s;  % endogenous
varexo ss;        % exogenous

% x : output gap (different from y: output)
% pi: rate of inflation
% i : nominal interest rate
% rr: real interest rate (Fisher equation)
% s : markp-shock
% ss: markup-shock innovation

parameters betta sig phin tetap alpha epsi omega kapa rhos; % parameters

betta = 0.99;
sig   = 1;
phin  = 1;
tetap = 0.75;
alpha = 0.33;
epsi  = 4;
omega = (1-alpha)/(1-alpha+alpha*epsi);
kapa  = ((1-tetap)*(1-tetap*betta)/tetap)*omega*(sig+(phin+alpha)/(1-alpha));

rhos  = 0.6;

%% 2) MODEL SETTING  - Write the equations of the model, expressing all at time t
%                      Note the log-linearized model (no steady state and linear specification)

model(linear);

x = x(+1) - 1/sig*(i - pi(+1));  % New-Keynesian IS Curve (NKIS)

pi = betta*pi(+1) + kapa*x + s;  % New-Keynesiam Phillips Curve (NKPC)

s = rhos*s(-1) + ss;             % mark-up shock

% monetary policy instrument
rr = i - pi(+1);                 % Fisher equation

end;

%% 3) No INITIAL VALUES since I am working with a log-linearized model

%% 4) CHECK STEADY STATE DIFFERENCES - Differences in computed steady state

%steady; % avoid to include it since, under commitment, there is no steady state
resid;

%% 5) SHOCKS 

shocks;
var ss = 0.01^2;
end;

%stoch_simul(irf=40);  % IRFs under no monetary intervention

%% 6) OPTIMAL DISCRETIONARY MONETARY POLICY

planner_objective pi^2 + kapa/epsi*x^2;

ramsey_policy(instruments=(i),irf=13,planner_discount=betta) x pi i rr s;

% This is what happens in the situation in which it is not implemented the Taylor rule, but the Fisher equation

%--------------------------------------------------------------------------
