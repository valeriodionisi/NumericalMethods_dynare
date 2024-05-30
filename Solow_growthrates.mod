% The Deterministic Solow (1956) model with growth rates changes
% Valerio Dionisi, v.dionisi@campus.unimib.it
% Department of Economics, Management and Statistics - University of Milano-Bicocca

% substitute 'true' with 'false' to run shock in n rather than in g
@#ifndef TFP_growth
    @#define TFP_growth=true  
@#endif

%--------------------------------------------------------------------------

%% 1) PREAMBLE - Parameter Calibration 
 
var c i y k exp_c exp_i exp_y exp_k             % endogenous variables: consumption, interest rate, output and capital
    g_k_aggregate g_k_per_capita g_k_intensive; % endogenous variables: growth rates

varexo A g n; %  exogenous variables: total factor productivity, productivity growth and population growth

% Parameters of the model, including variables at the steady state (ss)
parameters alpha delta s css iss yss kss Ass gss nss;

alpha = 0.33;  % capital share
delta = 0.025; % capital depreciation
s     = 0.2;   % savings rate

gss = 0.1;     % technology growth rate, steady state
nss = 0.2;     % population growth rate, steady state
Ass = 1;       % technology, steady state

css = (1 - s)*Ass*((Ass*s) / (delta + gss + nss))^(alpha / (1 - alpha));
iss = s*Ass*(Ass*s / (delta + gss + nss))^(alpha / (1 - alpha));
yss = Ass*(Ass*s / (delta + gss + nss))^(alpha / (1 - alpha));
kss = (Ass*s / (delta + gss + nss))^(1 / (1 - alpha));

%% 2) MODEL SETTING - Write the equations of the model, expressing all at time t (“stock at the end of the period” concept)
 
model;  % it declares the model
 
y = A*k(-1)^(alpha);               % production function. If the law of motion of capital
                                   % is at time t, then set k(-1) in the production function

k = i + (1 - delta - g - n)*k(-1); % law of motion for capital

y = c + i;                         % market clearing condition

c = (1 - s)*y;                     % consumption function

% exp values (to capture percentage changes)
exp_y = exp(y);                    % exp output
exp_c = exp(c);                    % exp consumption
exp_i = exp(i);                    % exp investment
exp_k = exp(k);                    % exp capital

% Growth rates
g_k_intensive  = exp(k)-exp(k(-1)); % capital growth rate between today and tomorrow
g_k_per_capita = g_k_intensive+g;   % capital per capita growth rate between today and tomorrow
g_k_aggregate  = g_k_intensive+g+n; % aggregate capital growth rate between today and tomorrow

end;  % end of model declaration

%% 3) MODEL AT THE STEADY STATE - Initialization of steady state values (copy the parameters an drop ss)

initval;  % initial values

% Recall that exogenous parameters have been already defined:
% alpha = 0.33;
% delta = 0.025;
% s = 0.2;

g = gss; % technology growth rate, initial value
n = nss; % population growth rate, initial value
A = Ass; % technology, initial value

c = css; % consumption, initial value
i = iss; % investment, initial value
y = yss; % output, initial value
k = kss; % capital, initial value

% otherwise, write
%c = (1 - s)*A*((A*s) / (delta + g + n))^(alpha / (1 - alpha));
%i = s*A*(A*s / (delta + g + n))^(alpha / (1 - alpha));
%y = A*(A*s / (delta + g + n))^(alpha / (1 - alpha));
%k = (A*s / (delta + g + n))^(1 / (1 - alpha));

% exp values (to capture percentage changes)
exp_y = exp(y);   % exp output
exp_c = exp(c);   % exp consumption
exp_i = exp(i);   % exp investment
exp_k = exp(k);   % exp predetermined capital

% Growth rates
g_k_intensive  = 0;
g_k_per_capita = g_k_intensive + g; 
g_k_aggregate  = g_k_intensive + g + n;

end;  % end of initial values declaration

%% 4) CHECK STEADY STATE DIFFERENCES - Differences in computed steady state

steady; % dynare computes the steady state. Steady state of the model for all the endogenous variables, assuming that exogenous variables are kept constant at the value declared in the initval block
resid;  % difference between model at the steady state and dynare steady state

% COMMENTS: no differences! in the second column, there are all zeros

%--------------------------------------------------------------------------

% Here, I am going to compute a permanent shock in g

@#if TFP_growth

%% 5) PERMANENT SHOCK ON TECHNOLOGY GROWTH RATE, g

initval;  % initial value declaration
g = 0.1;
end;
steady;   % steady states at the initial value

endval;   % final value declaration (i.e., the permanent value of the shock)
g = 0;
end;
steady;   % steady states at the final value

perfect_foresight_setup(periods = 100);
perfect_foresight_solver;

rplot A, c, i, y, k;

subplot(3,3,1)
plot(exp_c)
axis tight
ylabel('c')
title('PERMANENT SHOCK ON TECHNOLOGY GROWTH RATE, g')

subplot(3,3,2)
plot(exp_i)
axis tight
ylabel('i')

subplot(3,3,3)
plot(exp_y)
axis tight
ylabel('y')

subplot(3,3,4)
plot(exp_k)
axis tight
ylabel('k')

subplot(3,3,5)
plot(g_k_intensive)
axis tight
title('Growth rate of k in intensive form')
    
subplot(3,3,6)
plot(g_k_per_capita)
axis tight
title('Growth rate of per capita k')
    
subplot(3,3,7)
plot(g_k_aggregate)
axis tight
title('Growth rate of aggregate k')

%--------------------------------------------------------------------------

% Here, I am going to compute a permanent shock in n

@#else

%% 6) PERMANENT SHOCK ON POPULATION GROWTH RATE, n

initval;  % initial value declaration
n = 0.2;
end;
steady;   % steady states at the initial value

endval;   % final value declaration (i.e., the permanent value of the shock)
n = 0;
end;
steady;   % steady states at the final value

perfect_foresight_setup(periods = 100);
perfect_foresight_solver;

rplot A, c, i, y, k;

subplot(3,3,1)
plot(exp_c)
axis tight
ylabel('c')
title('PERMANENT SHOCK ON POPULATION GROWTH RATE, n')

subplot(3,3,2)
plot(exp_i)
axis tight
ylabel('i')

subplot(3,3,3)
plot(exp_y)
axis tight
ylabel('y')

subplot(3,3,4)
plot(exp_k)
axis tight
ylabel('k')

subplot(3,3,5)
plot(g_k_intensive)
axis tight
title('Growth rate of k in intensive form')
    
subplot(3,3,6)
plot(g_k_per_capita)
axis tight
title('Growth rate of per capita k')
    
subplot(3,3,7)
plot(g_k_aggregate)
axis tight
title('Growth rate of aggregate k')

%--------------------------------------------------------------------------

@#endif
