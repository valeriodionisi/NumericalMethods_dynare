% The Deterministic Ramsey-Cass-Koopmans model
% Valerio Dionisi, v.dionisi@campus.unimib.it
% Department of Economics, Management and Statistics - University of Milano-Bicocca

/*
 * This file uses Dynare's perfect foresight solver to study the transition 
 * behavior of simple non-stationary Ramsey-Cass-Koopmans economy with Cobb-Douglas production function
 * to its balanced growth path (BGP). It starts the simulation with a capital stock 
 * corresponding to 90% of its BGP value.
 * 
 * Notes:
 *  - The Ramsey-Cass-Koopmans model is solved here in aggregate, i.e. non-detrended form. Because
 *    in aggregate form there is no stationary steady state (only a BGP), one cannot
 *    use the steady-command after initval and endval to compute a conditional steady state.
 *  - The initial and terminal conditions are computed by using the analytical steady state
 *    in intensive form and then transform these values back to the aggregate BGP values
 *
 * This implementation was written by Johannes Pfeifer. 
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

@#define simulation_periods=100

%% 1) PREAMBLE - Parameter Calibration

var C I Y K g_K_aggregate g_K_per_capita g_K_intensive g_Y_aggregate g_Y_per_capita g_Y_intensive;

varexo A L; % Labor augmenting technology is growing in our experiment

parameters alpha beta delta g n;

alpha = 0.33; % capital share
delta = 0.1;  % capital depreciation
beta  = 0.99; % household discount factor
n     = 0.01; % population growth rate
g     = 0.02; % technology growth rate

%% 2) MODEL SETTING

model;

K = (1 - delta)*K(-1) + I;                        % law of motion for capital

I + C = Y;                                        % resource constraint

1/C = beta*1/C(+1)*(alpha*Y(+1)/K + (1 - delta)); % Euler (behavioral rule savings)

Y = K(-1)^alpha*(A*L)^(1 - alpha);                % production function

g_K_aggregate  = (K-K(-1))/K(-1);                                     % aggregate capital growth rate
g_K_per_capita = (K/L-K(-1)/L(-1))/(K(-1)/L(-1));                     % capital per capita growth rate
g_K_intensive  = (K/(A*L)-K(-1)/(A(-1)*L(-1)))/(K(-1)/(A(-1)*L(-1))); % capital growth rate
g_Y_aggregate  = (Y-Y(-1))/Y(-1);                                     % aggregate output growth rate
g_Y_per_capita = (Y/L-Y(-1)/L(-1))/(Y(-1)/L(-1));                     % output per capita growth rate
g_Y_intensive  = (Y/(A*L)-Y(-1)/(A(-1)*L(-1)))/(Y(-1)/(A(-1)*L(-1))); % output growth rate

end;

%% 3) MODEL AT THE STEADY STATE - Initialization of steady state values (copy the parameters an drop ss)

initval;

A = 1*(1+g)^0; % A_0 = 1
L = 1*(1+n)^0; % L_0 = 1

% starts the simulation with a capital stock corresponding to 50% of its BGP value
% steady state capital in intensive form multiplied by A*L
K = 0.5*(A*L)*(1+n)*(1+g)*((1/beta*(1+n)*(1+g)-(1-delta))/(alpha))^(1/(alpha-1));
Y = (K/(1 + n + g + n*g))^alpha*(A*L)^(1-alpha); %compute Y based on A_0, L_0, and K_(-1)
I = (1 - (1 - delta)/(1 + n + g + n*g))*K;
C = Y - I;
    
g_K_aggregate  = n + g + n*g;
g_K_per_capita = g; 
g_K_intensive  = 0; 
g_Y_aggregate  = n + g + n*g;
g_Y_per_capita = g; 
g_Y_intensive  = 0;

end;

%% 4) CHECK STEADY STATE DIFFERENCES - Differences in computed steady state

steady; % dynare computes the steady state. Steady state of the model for all the endogenous variables, assuming that exogenous variables are kept constant at the value declared in the initval block
resid;  % difference between model at the steady state and dynare steady state

%% 5) SHOCKS (define path of exogenous variables)

shock_vals_A = cumprod((1+g)*ones(@{simulation_periods},1))
shock_vals_L = cumprod((1+n)*ones(@{simulation_periods},1))

shocks;

var A;
periods 1:@{simulation_periods};
values (shock_vals_A);

var L;
periods 1:@{simulation_periods};
values (shock_vals_L);

end;

% set terminal condition to steady state value
endval;

A = 1*(1+g)^(@{simulation_periods}+1);
L = 1*(1+n)^(@{simulation_periods}+1);

% steady state capital in intensive form multiplied by A*L
K = (A*L)*(1 + n)*(1 + g)*((1/beta*(1 + n)*(1 + g) - (1 - delta))/(alpha))^(1/(alpha - 1));
Y = (K/(1 + n + g + n*g))^alpha*(A*L)^(1 - alpha);
I = (1 - (1 - delta)/(1 + n + g + n*g))*K;
C = Y - I;

g_K_aggregate  = n + g + n*g;
g_K_per_capita = g;
g_K_intensive  = 0;
g_Y_aggregate  = n + g + n*g;
g_Y_per_capita = g;
g_Y_intensive  = 0;

end;

perfect_foresight_setup(periods=@{simulation_periods});
perfect_foresight_solver; % compute the solution

% display simulation results
rplot K C Y;
rplot g_K_aggregate g_K_per_capita g_K_intensive;
rplot g_Y_aggregate g_Y_per_capita g_Y_intensive;
