function [K, P] = computeGainsMethodI(A, B, Q, R, P0)
% computeGainsMethodI - Method I: periodic gains by backward Riccati recursion
%
% Sets P_p = P_0 and recurses backward through the periodic Riccati equation
%
% Inputs:
%    A, B, Q, R - Periodic system and weights, cell arrays
%    P0         - Lifted Riccati solution
%
% Outputs:
%    K - Periodic gains, cell array
%    P - Periodic Riccati solution, cell array: P{l+1} = P_l

p = numel(A);
K = cell(1, p);
P = cell(1, p);
Pp1 = P0; % starts at P_p = P_0

% Retrieve solutions
for l = p:-1:1
    K{l} = (R{l} + B{l}' * Pp1 * B{l}) \ (B{l}' * Pp1 * A{l});
    P{l} = Q{l} + A{l}' * Pp1 * (A{l} - B{l} * K{l});
    Pp1 = P{l};
end

end
