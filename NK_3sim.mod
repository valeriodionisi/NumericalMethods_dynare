% The New Keynesian (NK) model with 3 equations
% Valerio Dionisi, v.dionisi@campus.unimib.it
% Department of Economics, Management and Statistics - University of Milano-Bicocca

%--------------------------------------------------------------------------

% 1) PREAMBLE - Parameter Calibration 

var y pi r s_D s_S s_r;                            % endogenous
varexo e_D e_S e_R;                                % exogenous

parameters beta sigma phi chi phi_pi phi_y lambda  % parameters
           rho_D rho_S rho_R kappa;

beta   = 0.99; 
sigma  = 1; 
phi    = 1; 
theta  = 1.0*3/4;
phi_pi = 1.5; 
phi_y  = 0.5/4; 
kappa  = ((1-theta)*(1-theta*beta)/theta)*(sigma+phi);

rho_D  = 0.9; 
rho_S  = 0.9; 
rho_R  = 0.3;

%% 2) MODEL SETTING  - Write the equations of the model, expressing all at time t
%                      Note the log-linearized model (no steady state and linear specification)

model(linear);  % it declares the model

y = y(+1) - 1/sigma*(r-pi(+1) - s_D); % New-Keynesian IS Curve (NKIS)

pi = beta*pi(+1) + kappa*y + s_S;     % New-Keynesiam Phillips Curve (NKPC)

r = phi_pi*pi + phi_y*y + s_r;        % Taylor rule

s_D = rho_D*s_D(-1) + e_D;            % demand (IS) shock
s_S = rho_S*s_S(-1) + e_S;            % supply (inflation) shock
s_r = rho_R*s_r(-1) + e_R;            % monetary policy shock

end;

%% 3) No INITIAL VALUES since I am working with a log-linearized model

%% 4) CHECK STEADY STATE DIFFERENCES - Differences in computed steady state
steady;
resid;

%% 5) SHOCKS 

shocks;

var e_D; stderr 1; 
var e_S; stderr 1; 
var e_R; stderr 1; 

end;

stoch_simul(irf=20) y r pi;

%--------------------------------------------------------------------------
