function [x, u, cost] = simulateLinear(A, B, Q, R, x0, gains, law, nSteps, sigma)
% simulateLinear - Closed loop on the linear periodic model, with either feedback law
%
% "stepwise": u_k = -K_l*x_k. "lifted": u_k = -Khat_l*x_s, with x_s the state at the start of the
% period.
%
% Inputs:
%    A, B, Q, R - Periodic system and weights, cell arrays
%    x0         - Initial state
%    gains      - K (cell) for "stepwise", Khat for "lifted"
%    law        - "stepwise" or "lifted"
%    nSteps     - Number of steps to simulate
%    sigma      - Standard deviation of the input noise (0 for a noise-free run)
%
% Outputs:
%    x    - States
%    u    - Commanded inputs
%    cost - Accumulated LQR cost

arguments
    A (1,:) cell
    B (1,:) cell
    Q (1,:) cell
    R (1,:) cell
    x0 (:,1) double
    gains
    law (1,1) string {mustBeMember(law, ["stepwise", "lifted"])}
    nSteps (1,1) {mustBeInteger, mustBePositive}
    sigma (1,1) double {mustBeNonnegative} = 0
end

p = numel(A);
m = size(B{1}, 2);
inputNoise = sigma * randn(m, nSteps);
x = [x0, zeros(numel(x0), nSteps)];
u = zeros(m, nSteps);
cost = 0;

for k = 1:nSteps
    l = mod(k - 1, p) + 1;
    if law == "stepwise"
        u(:, k) = -gains{l} * x(:, k);
    else
        if l == 1
            % Reset feedback state to start of period
            xStart = x(:, k);
        end
        u(:, k) = -gains((l - 1) * m + (1:m), :) * xStart;
    end
    cost = cost + x(:, k)' * Q{l} * x(:, k) + u(:, k)' * R{l} * u(:, k);
    x(:, k + 1) = A{l} * x(:, k) + B{l} * (u(:, k) + inputNoise(:, k));
end

end
