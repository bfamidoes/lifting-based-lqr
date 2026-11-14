function [K, Khat, P0] = computePeriodicLQR(A, B, Q, R, options)
% computePeriodicLQR - Periodic LQR gains computed from the lifted LQR solution
%
% Lifts the system over one period, solves the lifted Riccati equation, and computes the periodic
% gains by Method I (default, allows singular A_l) or Method II.
%
% Inputs:
%    A, B   - State and input matrices, cell arrays (cell l holds step l-1)
%    Q, R   - State and input weights, cell arrays
%    method - (name-value) 1 for Method I (default), 2 for Method II
%
% Outputs:
%    K    - Periodic gains, cell array: u_l = -K{l+1}*x_l
%    Khat - Nonzero block of the lifted gain
%    P0   - Lifted Riccati solution, equal to P_0 of the periodic Riccati equation

arguments
    A (1,:) cell
    B (1,:) cell
    Q (1,:) cell
    R (1,:) cell
    options.method (1,1) {mustBeMember(options.method, [1, 2])} = 1
end

[Ahat, Bhat, Qhat, Rhat, Shat] = liftPeriodicSystem(A, B, Q, R);
[P0, Khat] = solveLiftedRiccati(Ahat, Bhat, Qhat, Rhat, Shat);

if options.method == 1
    [K, ~] = computeGainsMethodI(A, B, Q, R, P0);
else
    K = computeGainsMethodII(A, B, Khat);
end

end
