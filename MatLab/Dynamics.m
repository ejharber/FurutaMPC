function [dxdt, y, A, B, C, D] = Dynamics(x, u)
% dxdt is the derivative of the states.
% [A B C D] are state space matrices linearized at the current operating point.

%% parameters
m_cart = 1;  % cart mass
m_pend = 1;  % pendulum mass
g = 9.81;    % gravity of earth
l = 1;       % pendulum length
r = 0.5;     %length of actuator

%% Obtain x, u and y
phi = x(1);
dphi = x(2);
theta = x(3);
dtheta = x(4);
% u
F = u;
% y
y = x;

%% Compute dxdt
dxdt = x;
% phi_dot
dxdt(1) = dphi;
% phi_dot_dot
dxdt(2) = (2*F*l + g*l*m_pend*r*sin(2*theta) - 2*dphi*dtheta*l^3*m_pend*sin(2*theta) - 2*dtheta^2*l^2*m_pend*r*sin(theta) + 2*dphi^2*m_pend*r^3*sin(theta)*(sin(theta)^2 - 1) + 2*dphi*dtheta*l*m_pend*r^2*sin(2*theta) - 2*dphi^2*l^2*m_pend*r*sin(theta)*(sin(theta)^2 - 1))/(2*l^3*m_pend + 2*l*m_cart*r^2 - 2*l^3*m_pend*cos(theta)^2);
% theta_dot
dxdt(3) = dtheta;
% theta_dot_dot
dxdt(4) = (2*F*l*r*cos(theta) + dphi^2*l^4*m_pend*sin(2*theta) - dphi^2*m_cart*r^4*sin(2*theta) + 2*g*l^3*m_pend*sin(theta) - 2*g*l^3*m_pend*(sin(theta) - sin(theta)^3) + 2*g*l*m_cart*r^2*sin(theta) + 2*g*l*m_pend*r^2*(sin(theta) - sin(theta)^3) + dphi^2*l^2*m_cart*r^2*sin(2*theta) - dphi^2*l^2*m_pend*r^2*sin(2*theta) - dtheta^2*l^2*m_pend*r^2*sin(2*theta) - 2*dphi^2*l^4*m_pend*cos(theta)^3*sin(theta) - 2*dphi^2*m_pend*r^4*cos(theta)^3*sin(theta) + 4*dphi^2*l^2*m_pend*r^2*cos(theta)^3*sin(theta) - 4*dphi*dtheta*l*m_pend*r^3*sin(theta)*(sin(theta)^2 - 1) + 4*dphi*dtheta*l^3*m_pend*r*sin(theta)*(sin(theta)^2 - 1))/(2*l^4*m_pend + 2*l^2*m_cart*r^2 - 2*l^4*m_pend*cos(theta)^2);

%% Obtain A/B/C/D from Jacobian
% Matricies derived in a seperate file
% LTI
A = [0, 1, 0, 0; ...
0, -(2*dtheta*m_pend*sin(2*theta)*l^3 - 2*dphi*m_pend*sin(2*theta)*cos(theta)*l^2*r - 2*dtheta*m_pend*sin(2*theta)*l*r^2 + 2*dphi*m_pend*sin(2*theta)*cos(theta)*r^3)/(2*l*(l^2*m_pend + m_cart*r^2 - l^2*m_pend*cos(theta)^2)), - (2*g*l*m_pend*r*sin(theta)^2 - 2*g*l*m_pend*r*cos(theta)^2 + 2*dphi^2*m_pend*r^3*cos(2*theta)*cos(theta) - dphi^2*m_pend*r^3*sin(2*theta)*sin(theta) + 4*dphi*dtheta*l^3*m_pend*cos(2*theta) + 2*dtheta^2*l^2*m_pend*r*cos(theta) - 4*dphi*dtheta*l*m_pend*r^2*cos(2*theta) - 2*dphi^2*l^2*m_pend*r*cos(2*theta)*cos(theta) + dphi^2*l^2*m_pend*r*sin(2*theta)*sin(theta))/(2*l*(l^2*m_pend + m_cart*r^2 - l^2*m_pend*cos(theta)^2)) - (l*m_pend*cos(theta)*sin(theta)*(m_pend*sin(2*theta)*cos(theta)*dphi^2*l^2*r - m_pend*sin(2*theta)*cos(theta)*dphi^2*r^3 - 2*m_pend*sin(2*theta)*dphi*dtheta*l^3 + 2*m_pend*sin(2*theta)*dphi*dtheta*l*r^2 - 2*m_pend*sin(theta)*dtheta^2*l^2*r + 2*g*m_pend*cos(theta)*sin(theta)*l*r + 2*F*l))/(l^2*m_pend + m_cart*r^2 - l^2*m_pend*cos(theta)^2)^2, -(2*dphi*m_pend*sin(2*theta)*l^3 + 4*dtheta*m_pend*sin(theta)*l^2*r - 2*dphi*m_pend*sin(2*theta)*l*r^2)/(2*l*(l^2*m_pend + m_cart*r^2 - l^2*m_pend*cos(theta)^2));...
0, 0, 0, 1;...
0, -(2*dphi*m_cart*r^4*sin(2*theta) - 2*dphi*l^4*m_pend*sin(2*theta) - 2*dphi*l^2*m_cart*r^2*sin(2*theta) + 2*dphi*l^2*m_pend*r^2*sin(2*theta) + 2*dphi*l^4*m_pend*sin(2*theta)*cos(theta)^2 + 2*dphi*m_pend*r^4*sin(2*theta)*cos(theta)^2 - 4*dphi*l^2*m_pend*r^2*sin(2*theta)*cos(theta)^2 - 2*dtheta*l*m_pend*r^3*sin(2*theta)*cos(theta) + 2*dtheta*l^3*m_pend*r*sin(2*theta)*cos(theta))/(2*l^2*(l^2*m_pend + m_cart*r^2 - l^2*m_pend*cos(theta)^2)), (2*dphi^2*l^4*m_pend*cos(2*theta) - 2*F*l*r*sin(theta) - 2*g*l^3*m_pend*cos(theta)^3 - 2*dphi^2*m_cart*r^4*cos(2*theta) + 2*g*l^3*m_pend*cos(theta) + 4*g*l^3*m_pend*cos(theta)*sin(theta)^2 + 2*g*l*m_cart*r^2*cos(theta) + 2*dphi^2*l^2*m_cart*r^2*cos(2*theta) - 2*dphi^2*l^2*m_pend*r^2*cos(2*theta) - 2*dtheta^2*l^2*m_pend*r^2*cos(theta)^2 + 2*dtheta^2*l^2*m_pend*r^2*sin(theta)^2 + 2*g*l*m_pend*r^2*cos(theta)^3 - 2*dphi^2*l^4*m_pend*cos(2*theta)*cos(theta)^2 - 2*dphi^2*m_pend*r^4*cos(2*theta)*cos(theta)^2 + 2*dphi^2*l^4*m_pend*sin(2*theta)*cos(theta)*sin(theta) + 2*dphi^2*m_pend*r^4*sin(2*theta)*cos(theta)*sin(theta) - 4*g*l*m_pend*r^2*cos(theta)*sin(theta)^2 + 4*dphi^2*l^2*m_pend*r^2*cos(2*theta)*cos(theta)^2 - 2*dphi*dtheta*l*m_pend*r^3*sin(2*theta)*sin(theta) + 2*dphi*dtheta*l^3*m_pend*r*sin(2*theta)*sin(theta) - 4*dphi^2*l^2*m_pend*r^2*sin(2*theta)*cos(theta)*sin(theta) + 4*dphi*dtheta*l*m_pend*r^3*cos(2*theta)*cos(theta) - 4*dphi*dtheta*l^3*m_pend*r*cos(2*theta)*cos(theta))/(2*l^2*(l^2*m_pend + m_cart*r^2 - l^2*m_pend*cos(theta)^2)) - (m_pend*cos(theta)*sin(theta)*(2*F*l*r*cos(theta) + dphi^2*l^4*m_pend*sin(2*theta) - dphi^2*m_cart*r^4*sin(2*theta) + 2*g*l^3*m_pend*sin(theta) - 2*g*l^3*m_pend*cos(theta)^2*sin(theta) + 2*g*l*m_cart*r^2*sin(theta) + dphi^2*l^2*m_cart*r^2*sin(2*theta) - dphi^2*l^2*m_pend*r^2*sin(2*theta) - dphi^2*l^4*m_pend*sin(2*theta)*cos(theta)^2 - dphi^2*m_pend*r^4*sin(2*theta)*cos(theta)^2 + 2*g*l*m_pend*r^2*cos(theta)^2*sin(theta) + 2*dphi^2*l^2*m_pend*r^2*sin(2*theta)*cos(theta)^2 - 2*dtheta^2*l^2*m_pend*r^2*cos(theta)*sin(theta) + 2*dphi*dtheta*l*m_pend*r^3*sin(2*theta)*cos(theta) - 2*dphi*dtheta*l^3*m_pend*r*sin(2*theta)*cos(theta)))/(l^2*m_pend + m_cart*r^2 - l^2*m_pend*cos(theta)^2)^2, -(2*dphi*m_pend*sin(2*theta)*cos(theta)*l^3*r + 4*dtheta*m_pend*cos(theta)*sin(theta)*l^2*r^2 - 2*dphi*m_pend*sin(2*theta)*cos(theta)*l*r^3)/(2*l^2*(l^2*m_pend + m_cart*r^2 - l^2*m_pend*cos(theta)^2))];
 
B = [0; 1/(l^2*m_pend + m_cart*r^2 - l^2*m_pend*cos(theta)^2); 0; (r*cos(theta))/(l*(l^2*m_pend + m_cart*r^2 - l^2*m_pend*cos(theta)^2))];
C = eye(4);
D = zeros(4,1);

