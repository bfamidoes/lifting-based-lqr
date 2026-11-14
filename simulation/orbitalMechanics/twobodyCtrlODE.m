function dydt = twobodyCtrlODE(~, y, mu, accelEci)
% twobodyCtrlODE - Two-body dynamics of chief and chaser, with a constant acceleration on the chaser
%
% The state stacks [r; v] of the chief, then [r; v] of the chaser, in ECI.
%
% Inputs:
%    ~        - Time (unused)
%    y        - State (12-by-1), m and m/s
%    mu       - Gravitational parameter (m^3/s^2)
%    accelEci - Chaser acceleration in ECI (m/s^2)
%
% Outputs:
%    dydt - State derivative (12-by-1)

dydt = zeros(12, 1);

% Chief
r1 = y(1:3);
v1 = y(4:6);
dydt(1:3) = v1;
dydt(4:6) = -mu * r1 / norm(r1)^3;

% Chaser
r2 = y(7:9);
v2 = y(10:12);
dydt(7:9)   = v2;
dydt(10:12) = -mu * r2 / norm(r2)^3 + accelEci;

end
