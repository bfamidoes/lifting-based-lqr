function [A, B, Q, R] = getAttitudeSystem(p)
% getAttitudeSystem - Euler-discretized magnetic attitude dynamics and LQR weights
%
% State [q_v; omega], input the magnetic dipole moment (A m^2); one orbit is divided into p steps.
%
% Inputs:
%    p - Steps per orbit
%
% Outputs:
%    A, B - State and input matrices, cell arrays
%    Q, R - State and input weights, cell arrays

J = [250; 150; 100];        % principal inertias, kg m^2
a = 7035e3;                 % orbit radius, m
mu = 3.986004418e14;        % Earth gravitational parameter, m^3/s^2
muF = 7.9e15;               % geomagnetic dipole strength, Wb m
iM = deg2rad(57);           % orbit inclination to the magnetic equator

w0 = sqrt(mu / a^3);
ts = 2 * pi / (w0 * p);     % Discretization time
b0 = muF / a^3;

Ac = [zeros(3), eye(3) / 2;
      8 * w0^2 * (J(3) - J(2)) / J(1), 0, 0, 0, 0, w0 * (- J(1) + J(2) - J(3)) / J(1);
      0, 6 * w0^2 * (J(3) - J(1)) / J(2), 0, 0, 0, 0;
      0, 0, 2 * w0^2 * (J(1) - J(2)) / J(3), w0 * (J(1) - J(2) + J(3)) / J(3), 0, 0];

A = repmat({eye(6) + Ac * ts}, 1, p);
Q = repmat({diag([2 2 2 1 1 1])}, 1, p);
R = repmat({1e-2 * eye(3)}, 1, p);
B = cell(1, p);
for l = 1:p
    theta = 2 * pi * (l - 1) / p;
    b = b0 * [sin(iM) * cos(theta); -cos(iM); 2 * sin(iM) * sin(theta)];
    B{l} = ts * [zeros(3); -[0, -b(3), b(2); b(3), 0, -b(1); -b(2), b(1), 0] ./ J];
end

end
