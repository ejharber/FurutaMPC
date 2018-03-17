clear;
clc;
clf;
close all;

% This furuta pendulum example is an extension of MATLAB's NMPC cart pole system
% Evan Harber
% eharber@andrew.cmu.edu

%% Product Requirement
% This example requires Optimization Toolbox(TM) to solve a nonlinear
% programming problem at each control interval.
if ~mpcchecktoolboxinstalled('optim')
    disp('Optimization Toolbox is required to run this example.')
    return
end

%% Control Objectives
% Assume the following initial conditions for the pendulum/cart assembly:
%
% * The actuator is stationary at _phi_ = |0|. 
%
% * The pendulum is in a downward equilibrium position where _theta_ = |-pi|.
%
% The control objectives are:
%
% * Swing-up control - Initially swing the pendulum up to an inverted
% equilibrium position where _phi_ = |0| and _theta_ = |0|.
%
% * Actuator position reference tracking - Move the cart to a new position with
% a step setpoint change, keeping the pendulum inverted.
%
% * Pendulum balancing - When an impulse disturbance of magnitude of |2|
% is applied to the inverted pendulum, keep the pendulum balanced and
% return the cart to its original position.

%% Setup
% Set the sample time:
Ts = 0.1;

% Set the prediction horizon: (How mant steps ahead the model looks)
N = 10;

x = [0;0;-pi;0];
uopt = zeros(N,1);
% x, the state space of the system follows - phi, dphi, theta, dtheta

% In the first stage of the simulation, the pendulum swings up from a
% downward equilibrium position to an inverted equilibrium position. The
% state references for this stage are all zero.
xref1 = [0;0;0;0];

% At a time of |10|, the cart moves from position |0| to |pi| and the
% references for the states become:
xref2 = [pi;0;0;0];

% In this example, to compute the optimal control sequence, use the
% |fmincon| function from the Optimization Toolbox as the nonlinear
% programming solver. Select the sequential quadratic programming (SQP)
% algorithm, which is one of the most common approaches in nonlinear MPC
% applications.
options = optimoptions('fmincon','Algorithm','sqp','Display','none');

% Run the simulation for |20| seconds, using |fmincon| with the defined
% nonlinear cost function and constraint function, |pendulumObjectiveFCN|
% and |pendulumConstraintFCN|, respectively.  The MV constraints are the
% lower and upper bounds of optimization variables in |fmincon|.
Duration = 20;

% Apply the MV constraints because the force has lower and upper bounds.
LB = -100*ones(N,1);
UB = 100*ones(N,1);

force = 0;
xHistory = x;

disp('Control Started');

tic 

t(1) = 0;
for ct = 1:(Duration/Ts)
    % Set references.
    if ct*Ts<10
        xref = xref1;
    else
        xref = xref2;
    end
    % Nonlinear MPC computation with full state feedback (no state estimator)
    COSTFUN = @(u) pendulumObjectiveFCN(u,x,Ts,N,xref,uopt(1));
    CONSFUN = @(u) pendulumConstraintFCN(u,x,Ts,N);
    uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,CONSFUN,options);    
    % Implement first optimal control move and update plant states.
    x = IntegrationEstimation(x, uopt(1), Ts, 30);
    % Save plant states for display.
    xHistory = [xHistory x];
    force(ct+1) = uopt(1);
    t(ct) = toc;
end
force 
disp('done')
figure(1);
subplot(2,2,1);
plot(0:Ts:Duration,xHistory(1,:));
xlabel('time (s)');
ylabel('phi');
title('actuator position');
subplot(2,2,2);
plot(0:Ts:Duration,xHistory(2,:));
xlabel('time (s)');
ylabel('dphi (1/s)');
title('actuator velocity');
subplot(2,2,3);
plot(0:Ts:Duration,xHistory(3,:));
xlabel('time (s)');
ylabel('theta');
title('pendulum angle');
subplot(2,2,4);
plot(0:Ts:Duration,xHistory(4,:));
xlabel('time (s)');
ylabel('thetadot (1/s)');
title('pendulum velocity');

figure(2);
plot(0:Ts:Duration,force);
xlabel('time (s)');
ylabel('force (N)');
title('Force vs. Time')

figure(3);
plot(t)
xlabel('iteration')
ylabel('time (s)')
title('Time vs. Iteration')
