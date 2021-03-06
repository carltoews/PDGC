clear all
close all

% Time interval and integration step-size
% (use course step-size for early investigation)
h = .01;
T = 10;
nSTEPS = T/h;

tspan = linspace(0, T, nSTEPS+1);

nControlPts = 51;

% Problem parameters
params.delta = .05;
params.p = 1;
params.q = .5;
params.c1 = .1;
params.c3 = 1;
params.r = @(t) 1 - .7*cos(2*pi*t);
params.k = 10;

ControlBounds = [-Inf Inf];
initState = [params.k; 4]; % The initial Effort value is free to be optimized by the solver
                           % The value of 4 serves as an initial guess only

% Construct Problem
prob = EngineeringProblem(params, ControlBounds);

% Use default PWLinearControl and RK4Integrator

% solve the problem
soln = single_shooting(prob, initState, tspan, nControlPts, ...
                       'FreeInitStates', 2, ... % Index of Effort state which has no fixed initial condition
                       'FreeStateBounds', [0 Inf]);

%% Plot the results
figure()
subplot(2,1,1);
plot(tspan, soln.x(tspan));
title('Optimal States');
legend('Fish population', 'Effort');

subplot(2,1,2);
plot(tspan, soln.u(tspan));
title('Optimal Control (dE/dt)');