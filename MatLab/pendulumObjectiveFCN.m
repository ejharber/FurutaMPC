function J = pendulumObjectiveFCN(u,x,Ts,N,xref,u0)
%% Cost function of nonlinear MPC for pendulum swing-up and balancing control
%
% Inputs:
%   u:      optimization variable, from time k to time k+N-1 
%   x:      current state at time k
%   Ts:     controller sample time
%   N:      prediction horizon
%   xref:   state references, constant from time k+1 to k+N
%   u0:     previous controller output at time k-1
%
% Output:
%   J:      objective function cost
%
% Copyright 2016 The MathWorks, Inc.

%% Nonlinear MPC design parameters
% Q matrix penalizes state deviations from the references.  Because we care
% more about cart position and pendulum angle, they have larger weights
% than their velocities.
Q = diag([10,1,10,1]);
% R matrix penalizes MV rate of change.  A small value is used here because
% aggressive control is desired.
R = 0.01;

%% Cost Calculation
% Set initial plant states, controller output and cost.
xk = x;
uk = u(1);
J = 0;
% Loop through each prediction step.
for ct=1:N
    % Obtain plant state at next prediction step.
    xk1 = IntegrationEstimation(xk, uk, Ts, 5);
    % accumulate state tracking cost from x(k+1) to x(k+N).
    J = J + (xk1-xref)'*Q*(xk1-xref);
    % accumulate MV rate of change cost from u(k) to u(k+N-1).
    if ct==1
        J = J + (uk-u0)'*R*(uk-u0);
    else
        J = J + (uk-u(ct-1))'*R*(uk-u(ct-1));
    end
    % Update xk and uk for the next prediction step.
    xk = xk1;
    if ct<N
        uk = u(ct+1);
    end
end

%% Cost Calculation for two inputs
% Set initial plant states, controller output and cost.
% xk = x;
% uk1 = u(1);
% uk2 = u(N+1);
% J = 0;
% % Loop through each prediction step.
% 
% for ct=1:N
%     % Obtain plant state at next prediction step.
%     xk1 = IntegrationEstimation(xk, uk, Ts, 5); %should use [uk1, uk2] as the input forces into the dynamic solver
%     % accumulate state tracking cost from x(k+1) to x(k+N).
%     J = J + pi-sqrt((xk1)'*Q*(xk1));
%     % accumulate MV rate of change cost from u(k) to u(k+N-1).
%     if ct==1
%         J = J + (uk1-u0(1))'*R*(uk1-u0(1)) + (uk2-u0(2))'*R*(uk2-u0(2));
%     else
%         J = J + (uk-u(ct-1))'*R*(uk-u(ct-1)) + ;
%     end
%     % Update xk and uk for the next prediction step.
%     xk = xk1;
%     if ct<N
%         uk = u(ct+1);
%     end
% end

%% Extend to 2 motors by having u be the size of 2*N w
