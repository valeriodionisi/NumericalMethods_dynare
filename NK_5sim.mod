% The New Keynesian (NK) model with 5 equations
% Valerio Dionisi, v.dionisi@campus.unimib.it
% Department of Economics, Management and Statistics - University of Milano-Bicocca

%--------------------------------------------------------------------------

%% 1) PREAMBLE - Parameter Calibration 

var y pi h i mu a nu;
varexo  eps_a eps_nu;

parameters beta theta phi sigma phi_pi;

beta   = 0.99;
theta  = 2/3;
phi    = 1;
sigma  = 1;
phi_pi = 1.75;

%% 2) MODEL SETTING  - Write the equations of the model, expressing all at time t
%                      Note the log-linearized model (no steady state and linear specification)

model(linear);

y = y(+1) - sigma^-1*(i-pi(+1));                      % New-Keynesian IS Curve (NKIS)

sigma*y + phi*h = a - mu;                             % labour supply

pi = beta*pi(+1) - (1-beta*theta)*(1-theta)/theta*mu; % New-Keynesiam Phillips Curve (NKPC)

y  = a + h;                                           % market clearing

i = phi_pi*pi + nu;                                   % Taylor rule

a = a(-1)*0.8 + eps_a;                                % mark-up shock
nu = nu(-1)*0.5 + eps_nu;                             % monetary policy shock

end;

%% 3) No INITIAL VALUES since I am working with a log-linearized model

%% 4) CHECK STEADY STATE DIFFERENCES - Differences in computed steady state
steady;
resid;

%% 5) SHOCKS 

shocks;

var eps_nu; stderr 1;
var eps_a; stderr 1;

end;

stoch_simul(irf=20);

%--------------------------------------------------------------------------
