% runAttitudeExample - Magnetic attitude control: lifted vs stepwise feedback under actuation noise
%
% Checks the optimal cost without noise, then compares both laws over ten orbits of actuation noise.

%% Setup

clear; clc;
setup;

p = 100;        % steps per orbit
x0 = [0.01; 0.01; 0.01; 1e-5; 1e-5; 1e-5]; % initial state

[A, B, Q, R] = getAttitudeSystem(p);
[K, Khat, P0] = computePeriodicLQR(A, B, Q, R);

%% Noise-free

% Without noise, both laws apply the same optimal inputs (with cost x0'*P0*x0)
nOrbitsNoiseFree = 50;
costOptimalNoiseFree = x0' * P0 * x0;

[~, ~, costStepwiseNoiseFree] = simulateLinear(A, B, Q, R, x0, K, "stepwise", nOrbitsNoiseFree * p);
[~, ~, costLiftedNoiseFree] = simulateLinear(A, B, Q, R, x0, Khat, "lifted", nOrbitsNoiseFree * p);
fprintf('Noise-free, %d orbits: x0''*P0*x0 = %.6e, stepwise = %.6e, lifted = %.6e\n', ...
    nOrbitsNoiseFree, costOptimalNoiseFree, costStepwiseNoiseFree, costLiftedNoiseFree);

%% Actuation noise

% Actuation noise on the commanded dipole; resetting the generator gives both laws the same realization
nOrbits = 10;
sigma = 0.025 / 3;
rng('default');
[xStepwise, uStepwise, costStepwise] = simulateLinear(A, B, Q, R, x0, K, "stepwise", nOrbits * p, sigma);
rng('default');
[xLifted, uLifted, costLifted] = simulateLinear(A, B, Q, R, x0, Khat, "lifted", nOrbits * p, sigma);
fprintf('Actuation noise, %d orbits: stepwise = %.4e, lifted = %.4e\n', nOrbits, costStepwise, costLifted);

%% Postprocess

plotAttitudeResults(xLifted, uLifted, xStepwise, uStepwise, p);