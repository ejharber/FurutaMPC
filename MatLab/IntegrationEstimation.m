function [xk1, yk] = IntegrationEstimation(xk, uk, Ts, numsteps)
%% Discrete-time nonlinear dynamic model of a pendulum on a cart at time k
%
% 4 states (xk): 
%   actuator position (phi)
%   actuator angular velocity (dphi): when positive, actuator moves to right
%   angle (theta): when 0, pendulum is at upright position
%   angular velocity (dtheta): when positive, pendulum moves anti-clockwisely
% 
% 1 inputs: (uk)
%   force (F): when positive, force pushes cart to right 
%
% 4 outputs: (yk)
%   same as states (i.e. all the states are measureable)
%
% xk1 is the states at time k+1.

%#codegen

% Repeat application of Euler method sampled at Ts/M.
M = numsteps;
delta = Ts/M;
xk1 = xk;
for ct=1:M
    if numsteps == 30
        k1 = Dynamics(xk1,uk);
        k2 = Dynamics(xk1 + k1*delta/2, uk);
        k3 = Dynamics(xk1 + k2*delta/2, uk);
        k4 = Dynamics(xk1 + k3*delta, uk);
        xk1 = xk1 + (k1/6 + k2/3 + k3/3 + k4/6)*delta;
    else
        xk1 = xk1 + delta*Dynamics(xk1,uk);
    end
end
yk = xk;
% Note that we choose the Euler method (first oder Runge-Kutta method)
% because it is more efficient for plant with non-stiff ODEs.  You can
% choose other ODE solvers such as ode23, ode45 for better accuracy or
% ode15s and ode23s for stiff ODEs.  Those solvers are available from
% MATLAB.
