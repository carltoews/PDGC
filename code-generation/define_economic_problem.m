%% Build symbolic functions
clear all
syms t delta p q E N J c1 c2 r k;
v = sym('v', [1 2]);
assume(v, 'real');

f(t, J, N, E) =  ...
   [r*N*(1 - N/k) - q*E*N; % State equation
   -exp(-delta*t)*(p*q*E*N - c1*E - c2*E^2)]; % Negation of Objective Equation (so solver can minimize)
   
v_times_f(t, J, N, E, v) = v*f;
dfdx_times_v = gradient(v_times_f, [N, J]);
dfdu_times_v = gradient(v_times_f, [E]);

%% Generate Code
matlabFunction(f, 'file', 'economic_F');
matlabFunction(dfdx_times_v, 'file', 'economic_dFdx_times_vec');
matlabFunction(dfdu_times_v, 'file', 'economic_dFdu_times_vec');