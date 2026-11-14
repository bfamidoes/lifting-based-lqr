function [Ac, Bc] = getThContinuousDynamics(f, e, mu, h)
% getThContinuousDynamics - Continuous Tschauner-Hempel matrices at a given true anomaly
%
% Derivatives are taken with respect to the true anomaly; units follow mu and h.
%
% Inputs:
%    f  - True anomaly (rad)
%    e  - Eccentricity
%    mu - Gravitational parameter
%    h  - Specific angular momentum
%
% Outputs:
%    Ac - State matrix (6-by-6), scaled TH coordinates
%    Bc - Input matrix (6-by-3), LVLH acceleration

kappa = 1 / (1 + e*cos(f));

Ac = [zeros(3), eye(3);
      3*kappa 0 0 0 2 0;
      0 0 0 -2 0 0;
      0 0 -1 0 0 0];

bCoeff = (h^(9/2)/mu^3) * kappa^3;
Bc = [zeros(3); bCoeff*eye(3)];

end
