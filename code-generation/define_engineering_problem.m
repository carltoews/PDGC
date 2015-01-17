%% Build symbolic functions
clear all
syms t delta p q E N J u c1 c3 r k;
v = sym('v', [1 3]);
assume(v, 'real');

f(t, N, E, J, u) =  ...
   [r*N*(1 - N/k) - q*E*N; % State equation
   u; % dEdt = u
   -exp(-delta*t)*(p*q*E*N - c1*E - c3*u^2)]; % J = -objective (so solver can minimize)
   
v_times_f(t, N, E, J, u, v) = v*f;
dfdx_times_v = gradient(v_times_f, [N, E, J]);
dfdu_times_v = gradient(v_times_f, [u]);

%% Generate Code
matlabFunction(f, 'file', 'engineering_F');
matlabFunction(dfdx_times_v, 'file', 'engineering_dFdx_times_vec');
matlabFunction(dfdu_times_v, 'file', 'engineering_dFdu_times_vec');