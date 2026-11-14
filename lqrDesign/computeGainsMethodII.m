function K = computeGainsMethodII(A, B, Khat)
% computeGainsMethodII - Method II: periodic gains extracted from the lifted gain
%
% Requires A_0, ..., A_(p-2) nonsingular.
%
% Inputs:
%    A, B - Periodic system, cell arrays
%    Khat - Nonzero block of the lifted gain
%
% Outputs:
%    K - Periodic gains, cell array

[n, m] = size(B{1});
p = numel(A);
K = cell(1, p);
PiProd = eye(n);
RCOND_TOL = 1e-12;
minRcond = 1;

% Retrieve gains
for l = 1:p
    KhatL = Khat((l - 1) * m + (1:m), :);
    minRcond = min(minRcond, rcond(PiProd));
    K{l} = KhatL / PiProd;
    PiProd = A{l} * PiProd - B{l} * KhatL;
end

if minRcond < RCOND_TOL
    warning('computeGainsMethodII:singularPi', ...
        'Product of Pi_l nearly singular (min rcond %.1e): the gains may be inaccurate; use Method I.', minRcond);
end

end
