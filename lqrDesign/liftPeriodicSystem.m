function [Ahat, Bhat, Qhat, Rhat, Shat] = liftPeriodicSystem(A, B, Q, R)
% liftPeriodicSystem - Lifted system and cost over one period, from the block rows of the lifted system
%
% Computes the matrices of the reduced lifted Riccati equation, whose solution is P_0.
%
% Inputs:
%    A, B, Q, R - Periodic system and weights, cell arrays
%
% Outputs:
%    Ahat, Bhat       - Lifted state and input matrices
%    Qhat, Rhat, Shat - Lifted weights
%
% Reference: Yang (2018), Automatica 87:383-388, doi:10.1016/j.automatica.2017.10.019

[n, m] = size(B{1});
p = numel(A);

% Allocation
Abar1 = zeros((p - 1) * n, n);
Bbar1 = zeros((p - 1) * n, p * m);
QAbar1 = zeros((p - 1) * n, n);         % Qbar1*Abar1
QBbar1 = zeros((p - 1) * n, p * m);     % Qbar1*Bbar1
Ahat = A{1};
Bhat = zeros(n, p * m);
Bhat(:, 1:m) = B{1};

% Block rows, propagated from x_1 to x_p
for l = 1:p - 1
    rows = (l - 1) * n + (1:n);
    Abar1(rows, :) = Ahat;
    Bbar1(rows, :) = Bhat;
    QAbar1(rows, :) = Q{l + 1} * Ahat;
    QBbar1(rows, :) = Q{l + 1} * Bhat;
    Ahat = A{l + 1} * Ahat;
    Bhat = A{l + 1} * Bhat;
    Bhat(:, l * m + (1:m)) = B{l + 1};
end

% Lifted weights
Qhat = Q{1} + Abar1' * QAbar1;
Rhat = blkdiag(R{:}) + Bbar1' * QBbar1;
Shat = Abar1' * QBbar1;

end
