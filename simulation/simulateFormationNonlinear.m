function [xLvlh, accel, cost] = simulateFormationNonlinear(chief0, chaser0, orbit, gains, law, Q, R, nOrbits)
% simulateFormationNonlinear - Formation feedback on the nonlinear two-body dynamics
%
% Each step converts the ECI states to the TH state, applies either feedback law (as in
% simulateLinear), and holds the commanded acceleration over one true-anomaly step.
%
% Inputs:
%    chief0, chaser0  - Initial [r; v] in ECI (m, m/s)
%    orbit            - Orbit data
%    gains, law       - As in simulateLinear
%    Q, R             - Weights, cell arrays
%    nOrbits          - Number of orbits
%
% Outputs:
%    xLvlh - Chaser state relative to the chief in LVLH (km, km/s)
%    accel - Commanded LVLH acceleration (km/s^2)
%    cost  - Accumulated LQR cost in TH coordinates

arguments
    chief0 (6,1) double
    chaser0 (6,1) double
    orbit struct
    gains
    law (1,1) string {mustBeMember(law, ["stepwise", "lifted"])}
    Q (1,:) cell
    R (1,:) cell
    nOrbits (1,1) {mustBePositive, mustBeInteger}
end

muSI = orbit.mu * 1e9;
p = numel(Q);
m = 3;
nSteps = nOrbits * p;
% Get time grid for anomaly steps of simulation
tGrid = arrayfun(@(f) trueAnomaly2Time(f, orbit.e, orbit.n), orbit.f0 + (0:nSteps) * orbit.df);
odeOptions = odeset('RelTol', 1e-13, 'AbsTol', 1e-15);

y = [chief0; chaser0];
xLvlh = zeros(6, nSteps + 1);
accel = zeros(m, nSteps);
cost = 0;

for k = 1:nSteps
    l = mod(k - 1, p) + 1;
    [pRel, vRel, eciToLvlhDcm] = eci2Lvlh(y(1:3), y(4:6), y(7:9), y(10:12));
    xLvlh(:, k) = [pRel; vRel] / 1e3;
    f = eci2TrueAnomaly(y(1:3), y(4:6), muSI);
    xTH = buildLvlh2Th(f, orbit.e, orbit.mu, orbit.h) * xLvlh(:, k);

    if law == "stepwise"
        u = -gains{l} * xTH;
    else
        if l == 1
            xStart = xTH;
        end
        u = -gains((l - 1) * m + (1:m), :) * xStart;
    end
    cost = cost + xTH' * Q{l} * xTH + u' * R{l} * u;

    accel(:, k) = u * orbit.uRef;
    accelEci = eciToLvlhDcm' * accel(:, k) * 1e3;
    % Propagate until feedback is refreshed
    [~, ySol] = ode45(@(t, yy) twobodyCtrlODE(t, yy, muSI, accelEci), tGrid(k:k+1), y, odeOptions);
    y = ySol(end, :)';
end

[pRel, vRel] = eci2Lvlh(y(1:3), y(4:6), y(7:9), y(10:12));
xLvlh(:, end) = [pRel; vRel] / 1e3;

end
