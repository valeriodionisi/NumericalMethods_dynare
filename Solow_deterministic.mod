% The Deterministic Solow (1956) model
% Valerio Dionisi, v.dionisi@campus.unimib.it
% Department of Economics, Management and Statistics - University of Milano-Bicocca



% substitute 'true' with 'false' to run shocks in g, n, and TFP shock in period 5
@#ifndef TFP_transitory_permanent
    @#define TFP_transitory_permanent=false  
@#endif

%--------------------------------------------------------------------------

% First, I simulate a TFP shock through A. Therefore, I should have only
% one exogenous variable, taking the rest as parameters.

@#if TFP_transitory_permanent

%% 1.1) PREAMBLE - Parameter Calibration 
 
var c i y k exp_c exp_i exp_y exp_k;  % endogenous variables: consumption, interest rate, output and capital

varexo A; %  exogenous variables: total factor productivity

% Parameters of the model, including variables at the steady state (ss)
parameters alpha delta s g n css iss yss kss Ass gss nss;

alpha = 0.33;  % capital share
delta = 0.025; % capital depreciation
s     = 0.2;   % savings rate
g     = 0.1;   % technology growth rate
n     = 0.2;   % population growth rate

Ass = 1;       % technology, steady state

css = (1 - s)*Ass*((Ass*s) / (delta + gss + nss))^(alpha / (1 - alpha));
iss = s*Ass*(Ass*s / (delta + gss + nss))^(alpha / (1 - alpha));
yss = Ass*(Ass*s / (delta + gss + nss))^(alpha / (1 - alpha));
kss = (Ass*s / (delta + gss + nss))^(1 / (1 - alpha));

%% 2.1) MODEL SETTING - Write the equations of the model, expressing all at time t (“stock at the end of the period” concept)
 
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
exp_k = exp(k(+1));                % exp capital

end;  % end of model declaration

%% 3.1) MODEL AT THE STEADY STATE - Initialization of steady state values (copy the parameters an drop ss)

initval;  % initial values

% Recall that exogenous parameters have been already defined:
% alpha = 0.33;
% delta = 0.025;
% s     = 0.2;
% g     = 0.1;
% n     = 0.2;

A = Ass;   % technology, initial value

c = css;   % consumption, initial value
i = iss;   % investment, initial value
y = yss;   % output, initial value
k = kss;   % capital, initial value

% otherwise, write steady_state_model; and thereafter
% c = (1 - s)*A*((A*s) / (delta + g + n))^(alpha / (1 - alpha));
% i = s*A*(A*s / (delta + g + n))^(alpha / (1 - alpha));
% y = A*(A*s / (delta + g + n))^(alpha / (1 - alpha));
% k = (A*s / (delta + g + n))^(1 / (1 - alpha));

% exp values (to capture percentage changes)
exp_y = exp(y);   % exp output
exp_c = exp(c);   % exp consumption
exp_i = exp(i);   % exp investment
exp_k = exp(k);   % exp predetermined capital

end;  % end of initial values declaration

%% 4.1) CHECK STEADY STATE DIFFERENCES - Differences in computed steady state

steady; % dynare computes the steady state. Steady state of the model for all the endogenous variables, assuming that exogenous variables are kept constant at the value declared in the initval block
resid;  % difference between model at the steady state and dynare steady state

% COMMENTS: no differences! in the second column, there are all zeros

%% 5) TRANSITORY SHOCK ON TFP

shocks;       % it identifies transitory shocks! "overwrite" replaces all the previous shocks blocks
var A;
periods 1:4;  % one year shock (remember: quarterly data!)
values 1.01;  %  0.01 shock in A. Remember that initial value is A = 1
end;

perfect_foresight_setup(periods = 25);  % number of periods for simulation
perfect_foresight_solver;

rplot A, c, i, y, k;

subplot(4,1,1)
plot(exp_c)
axis tight
ylabel('c')
title('ONE YEAR TRANSITORY SHOCK ON TFP')

subplot(4,1,2)
plot(exp_i)
axis tight
ylabel('i')

subplot(4,1,3)
plot(exp_y)
axis tight
ylabel('y')

subplot(4,1,4)
plot(exp_k)
axis tight
ylabel('k')

% COMMENTS: Fluctuations and then stabilization (See 00_.endo_simul)

%% 6) PERMANENT SHOCK ON TFP

% Under deterministic model, one should write two steady states: the first
% one is the initial (that with the calibrated value at the beginning); the
% second is the steady state after the shock occurs, that is, the steady
% state under the new value of the variable that experiences the shock.

initval;  % initial value declaration
A = 1;
end;
steady;   % steady states at the initial value

endval;   % final value declaration (i.e., the permanent value of the shock)
A = 1.01;
end;
steady;   % steady states at the final value

perfect_foresight_setup(periods = 100);
perfect_foresight_solver;

rplot A, c, i, y, k;

subplot(4,1,1)
plot(exp_c)
axis tight
ylabel('c')
title('PERMANENT SHOCK ON TFP')

subplot(4,1,2)
plot(exp_i)
axis tight
ylabel('i')

subplot(4,1,3)
plot(exp_y)
axis tight
ylabel('y')

subplot(4,1,4)
plot(exp_k)
axis tight
ylabel('k')

% COMMENTS: Everything increases permanently even by changing time span of
            the IRFs (See the command window)

%--------------------------------------------------------------------------

% Now, I am going to add shocks in g and n. The codes are the same as 
% before, with the only difference now these should appear as exogenous
% variables, and no more as parameters (adding their staedy-state values).

@#else

%% 1.2) PREAMBLE - Parameter Calibration 
 
var c i y k exp_c exp_i exp_y exp_k;  % endogenous variables: consumption, interest rate, output and capital

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

%% 2.2) MODEL SETTING - Write the equations of the model, expressing all at time t (“stock at the end of the period” concept)
 
model;  % it declares the model
 
y = A*k(-1)^(alpha);               % production function. If the law of motion of capital
                                   % is at time t, then set k(-1) in the production function

k = i + (1 - delta - g - n)*k(-1); % law of motion for capital

y = c + i;                         % market clearing condition

c = (1 - s)*y;                     % consumption function

% Log values
exp_y = exp(y);                    % exp output
exp_c = exp(c);                    % exp consumption
exp_i = exp(i);                    % exp investment
exp_k = exp(k(+1));                % exp capital

end;  % end of model declaration

%% 3.2) MODEL AT THE STEADY STATE - Initialization of steady state values (copy the parameters an drop ss)

initval;  % initial values

% Recall that exogenous parameters have been already defined:
% alpha = 0.33;
% delta = 0.025;
% s     = 0.2;

g = gss;   % technology growth rate, initial value
n = nss;   % population growth rate, initial value
A = Ass;   % technology, initial value

c = css;   % consumption, initial value
i = iss;   % investment, initial value
y = yss;   % output, initial value
k = kss;   % capital, initial value

% otherwise, write steady_state_model; and thereafter
% c = (1 - s)*A*((A*s) / (delta + g + n))^(alpha / (1 - alpha));
% i = s*A*(A*s / (delta + g + n))^(alpha / (1 - alpha));
% y = A*(A*s / (delta + g + n))^(alpha / (1 - alpha));
% k = (A*s / (delta + g + n))^(1 / (1 - alpha));

% Log values
exp_y = exp(y);                % exp output
exp_c = exp(c);                % exp consumption
exp_i = exp(i);                % exp investment
exp_k = exp(k);                % exp predetermined capital

end;  % end of initial values declaration

%% 4.2) CHECK STEADY STATE DIFFERENCES - Differences in computed steady state

steady; % dynare computes the steady state. Steady state of the model for all the endogenous variables, assuming that exogenous variables are kept constant at the value declared in the initval block
resid;  % difference between model at the steady state and dynare steady state

% COMMENTS: no differences! in the second column, there are all zeros

%% 7) PERMANENT CHANGE TO TFP OCCURRING IN PERIOD 5

initval;  % initial value declaration
A = 1;
end;
steady;   % steady states at the initial value

endval;   % final value declaration (i.e., the permanent value of the shock)
A = 1.01;
end;
steady;   % steady states at the final value

shocks;
var A;
periods 1:4;
values 1;
end;

perfect_foresight_setup(periods = 25);
perfect_foresight_solver;

rplot A, c, i, y, k;

subplot(4,1,1)
plot(c)
axis tight
ylabel('c')
title('PERMANENT CHANGE TO TFP IN PERIOD 5')

subplot(4,1,2)
plot(i)
axis tight
ylabel('i')

subplot(4,1,3)
plot(y)
axis tight
ylabel('y')

subplot(4,1,4)
plot(k)
axis tight
ylabel('k')

%% 8) PERMANENT CHANGES TO g AND n

% Initial values are g = 1.1, n = 1.2, and s = 0.2;

% Let's start from g ...
initval;  % initial value declaration
g = 0.1;
end;
steady;   % steady states at the initial value

endval;   % final value declaration (i.e., the permanent value of the shock)
g = 0.2;
end;
steady;   % steady states at the final value

perfect_foresight_setup(periods = 25);
perfect_foresight_solver;

rplot A, c, i, y, k;

subplot(4,1,1)
plot(c)
axis tight
ylabel('c')
title('PERMANENT CHANGE IN g')

subplot(4,1,2)
plot(i)
axis tight
ylabel('i')

subplot(4,1,3)
plot(y)
axis tight
ylabel('y')

subplot(4,1,4)
plot(k)
axis tight
ylabel('k')

% COMMENTS: Everything decreases (c = 0.269846; i = 0.0674616; y = 0.337308; k = 0.0278192)

% ... and then let's look at n
initval;  % initial value declaration
n = 0.2;
end;
steady;   % steady states at the initial value

endval;   % final value declaration (i.e., the permanent value of the shock)
n = 0.3;
end;
steady;   % steady states at the final value

perfect_foresight_setup(periods = 25);
perfect_foresight_solver;

rplot A, c, i, y, k;

subplot(4,1,1)
plot(c)
axis tight
ylabel('c')
title('PERMANENT CHANGE IN n')

subplot(4,1,2)
plot(i)
axis tight
ylabel('i')

subplot(4,1,3)
plot(y)
axis tight
ylabel('y')

subplot(4,1,4)
plot(k)
axis tight
ylabel('k')

% COMMENTS: Everything decreases (c = 0.264529; i = 0.0661321; y = 0.330661; k = 0.0261909)

%--------------------------------------------------------------------------

@#endif
