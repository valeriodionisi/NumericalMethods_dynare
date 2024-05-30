% The Real Business Cycle (RBC) model
% Valerio Dionisi, v.dionisi@campus.unimib.it
% Department of Economics, Management and Statistics - University of Milano-Bicocca

%--------------------------------------------------------------------------

% 1) PREAMBLE - Parameter Calibration 

var C K L w r A;
varexo e;

parameters alpha delta gamma rho rho_A g;

alpha = 0.33;  % capital share of output
delta = 0.1;   % capital depreciation
gamma = 0;     % inverse of Frisch elasticity
rho   = 0.03;  % discount rate
rho_A = 0.97;  % persistence of technology AR(1) process
g     = 0.015; % labor-augmenting productivity growth rate

% 2) MODEL SETTING  - Write the equations of the model, expressing all at time t (“stock at the end of the period” concept)

model;

1/C = 1/(1+rho)*(1/(C(+1)*(1+g)))*(r(+1)+1-delta);               % Euler

L^gamma = w/C;                                                   % labour supply

r = alpha*A*(K(-1)/(1+g))^(alpha-1)*L^(1-alpha);                 % interest rate

w = (1-alpha)*A*(K(-1)/(1+g))^alpha*L^(-alpha);                  % wage rate

K+C = (K(-1)/(1+g))*(1-delta)+A*(K(-1)/(1+g))^alpha*L^(1-alpha); % resource constraint

log(A) = rho_A*log(A(-1)) + e;                                   % technology dynamics

end;

%% 3) MODEL AT THE STEADY STATE  - Initialization of steady state values

steady_state_model;

A = 1;                                                % technology

r = (1+g)*(1+rho)+delta-1;                            % rental rate

L = ((1-alpha)/(r/alpha-delta-g))*r/alpha;            % labour

K = (1+g)*(r/alpha)^(1/(alpha-1))*L;                  % capital

C = (1-delta)*K/(1+g)+(K/(1+g))^alpha*L^(1-alpha)-K;  % resource constraint

w = C;                                                % wage rate

end;

%% 4) CHECK STEADY STATE DIFFERENCES - Differences in computed steady state
steady;
resid;

%% 5) SHOCK ON TFP

shocks;
var e; stderr 0.01; % technology innovation in AR(1)
end;

stoch_simul(order=2,irf=200);

% Note, if one includes "order" in stoch_simul(), it stands for the order of
% the Taylor expansion. The default order is 2.

%--------------------------------------------------------------------------
