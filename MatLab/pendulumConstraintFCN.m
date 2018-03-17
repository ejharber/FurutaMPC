function [c, ceq] = pendulumConstraintFCN(u,x,Ts,N)
% Inputs:
%   u:      optimization variable, from time k to time k+N-1 
%   x:      current state at time k
%   Ts:     controller sample time
%   N:      prediction horizon
%
% Output:
%   c:      inequality constraints applied across prediction horizon
%   ceq:    equality constraints (empty)

%% Nonlinear MPC design parameters
% Range of cart position: from -10 to 10
phiMin = -10;
phiMax = 10;

%% Inequality constraints calculation
c = zeros(N*2,1);
% Apply 2*N cart position constraints across prediction horizon, from time
% k+1 to k+N
xk = x;
uk = u(1);
for ct=1:N
    % obtain new cart position at next prediction step
    xk1 = IntegrationEstimation(xk, uk, Ts, 5);
    % -z + zMin < 0
    c(2*ct-1) = -xk1(1)+phiMin;
    % z - zMax < 0
    c(2*ct) = xk1(1)-phiMax;
    % update plant state and input for next step
    xk = xk1;
    if ct<N
        uk = u(ct+1);
    end
end
%% No equality constraints
ceq = [];

