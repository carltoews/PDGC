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
params.c2 = .1;
params.r = 1;
params.k = 10;

ControlBounds = [0 Inf];
initState = params.k;

% Construct Problem
prob = EconomicProblem(params, ControlBounds);

% Use default PWLinearControl and RK4Integrator

% solve the problem
soln = single_shooting(prob, initState, tspan, nControlPts);

% Plot the results
figure()
subplot(2,1,1);
plot(tspan, soln.x(tspan));
title('Fish population');

subplot(2,1,2);
plot(tspan, soln.u(tspan));
title('Optimal Effort');