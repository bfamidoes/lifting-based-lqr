function f = mean2TrueAnomaly(M, e)
% mean2TrueAnomaly - True anomaly from mean anomaly
%
% Solves Kepler's equation M = E - e*sin(E) for the eccentric anomaly E by Newton's method.
%
% Inputs:
%    M - Mean anomaly (rad)
%    e - Eccentricity (elliptic orbit)
%
% Outputs:
%    f - True anomaly (rad)

E = M + e * sin(M);
MAX_ITER = 1e3;
TOL = 1e-14;

% Simple Newton-Raphson iteration method
for iter = 1:MAX_ITER
    dE = (E - e * sin(E) - M) / (1 - e * cos(E));
    E = E - dE;
    if abs(dE) < TOL
        break
    end
end

if ~(abs(dE) < TOL)
    warning('mean2TrueAnomaly:noConvergence', ...
        'Kepler''s equation not solved within %d iterations (last step %.2e rad).', MAX_ITER, abs(dE));
end

f = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));

end
