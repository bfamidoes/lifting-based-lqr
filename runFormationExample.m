% runFormationExample - Formation flying on a Tundra orbit: lifted vs stepwise feedback
%
% Designs on the Euler-discretized Tschauner-Hempel model, then simulates both laws on the
% nonlinear two-body dynamics.

%% Setup

clear; clc;
setup;

p = 100;            % steps per orbit, uniform in true anomaly
nOrbits = 5;
uRef = 1e-7;        % input scaling c_N, km/s^2

% Keplerian elements [a (m), e, i, argument of perigee, RAAN, mean anomaly (deg)]
chiefElements = [42164e3, 0.27, 63.4, 270, 190, 0];
chaserElements = chiefElements + [-300, -2e-4, 4e-3, 0, 6e-3, 5e-3];

f0 = mean2TrueAnomaly(deg2rad(chiefElements(6)), chiefElements(2));
[A, B, orbit] = getTschaunerHempelSystem(p, chiefElements(1), chiefElements(2), f0, uRef);
Q = repmat({eye(6)}, 1, p);
R = repmat({diag([30, 30, 3])}, 1, p);
[K, Khat] = computePeriodicLQR(A, B, Q, R);

%% Simulation

chief0 = computeEciState(chiefElements, orbit.mu * 1e9);
chaser0 = computeEciState(chaserElements, orbit.mu * 1e9);
[xStepwise, ~, costStepwise] = simulateFormationNonlinear(chief0, chaser0, orbit, K, "stepwise", Q, R, nOrbits);
[xLifted, ~, costLifted] = simulateFormationNonlinear(chief0, chaser0, orbit, Khat, "lifted", Q, R, nOrbits);
fprintf('Nonlinear dynamics, %d orbits: stepwise = %.4f, lifted = %.4f, ratio lifted/stepwise = %.4f\n', ...
    nOrbits, costStepwise, costLifted, costLifted / costStepwise);

%% Postprocess

plotFormationResults(xLifted, xStepwise, p);

%% Local

function state = computeEciState(elements, mu)
% computeEciState - ECI state from Keplerian elements
%
% Retrieves ECI state.
%
% Inputs:
%    elements - [a (m), e, i, argument of perigee, RAAN, mean anomaly (deg)]
%    mu       - Gravitational parameter (m^3/s^2)
%
% Outputs:
%    state - [r; v] in ECI (m, m/s)

angles = deg2rad(elements(3:6));
f = mean2TrueAnomaly(angles(4), elements(2));
[r, v] = kepler2Eci(elements(1), elements(2), angles(1), angles(2), angles(3), f, mu);
state = [r; v];

end
