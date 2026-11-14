function [P0, Khat] = solveLiftedRiccati(Ahat, Bhat, Qhat, Rhat, Shat)
% solveLiftedRiccati - Lifted Riccati solution and lifted gain
%
% idare is used for readability and robustness; a direct QZ solution is faster.
%
% Inputs:
%    Ahat, Bhat, Qhat, Rhat, Shat - Lifted system and weights, from liftPeriodicSystem
%
% Outputs:
%    P0   - Stabilizing solution
%    Khat - Nonzero block of the lifted gain
%
% Reference: Yang (2018), Automatica 87:383-388, doi:10.1016/j.automatica.2017.10.019

% Note: A direct QZ solution is significantly faster for large periods. For
% this demo, with p = 100, idare is simpler to read and shows good performance
[P0, Khat] = idare(Ahat, Bhat, Qhat, Rhat, Shat, []);
if isempty(P0)
    error('solveLiftedRiccati:noSolution', 'idare found no stabilizing solution of the lifted Riccati equation');
end

end
