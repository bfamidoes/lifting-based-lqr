function plotAttitudeResults(xLifted, uLifted, xStepwise, uStepwise, p)
% plotAttitudeResults - Attitude, body rate and commanded dipole for both feedback laws
%
% One state figure per law (quaternion on top, body rate below), then one figure of the commanded
% dipole with the lifted law on top and the stepwise law below.
%
% Inputs:
%    xLifted, uLifted     - States and inputs under the lifted law, from simulateLinear
%    xStepwise, uStepwise - States and inputs under the stepwise law
%    p                    - Steps per orbit
%
% Outputs:
%    none (three figures)

nSteps = size(uLifted, 2);
orbit = (1:nSteps) / p;
laws = ["Lifted feedback", "Stepwise feedback"];
states = {xLifted(:, 2:end), xStepwise(:, 2:end)};
inputs = {uLifted, uStepwise};

% State plot for Lifted/Stepwise
for i = 1:2
    figure('Name', laws(i));
    colororder([0 0 1; 1 0 0; 0 1 0]);
    subplot(2, 1, 1);
    plot(orbit, states{i}(1:3, :));
    grid on;
    title(laws(i), 'Interpreter', 'latex');
    ylabel('Quaternion component', 'Interpreter', 'latex');
    legend('$q_1$', '$q_2$', '$q_3$', 'Interpreter', 'latex', 'Location', 'best');
    subplot(2, 1, 2);
    plot(orbit, states{i}(4:6, :));
    grid on;
    ylabel('Body rate (rad/s)', 'Interpreter', 'latex');
    legend('$\omega_1$', '$\omega_2$', '$\omega_3$', 'Interpreter', 'latex', 'Location', 'best');
    xlabel('Orbit', 'Interpreter', 'latex');
end

% Control plot
figure('Name', 'Commanded dipole');
colororder([0 0 1; 1 0 0; 0 1 0]);
for i = 1:2
    subplot(2, 1, i);
    plot(orbit, inputs{i});
    grid on;
    title(laws(i), 'Interpreter', 'latex');
    ylabel({'Magnetic moment', '(A$\cdot$m$^2$)'}, 'Interpreter', 'latex');
    legend('$m_1$', '$m_2$', '$m_3$', 'Interpreter', 'latex', 'Location', 'best');
end
xlabel('Orbit', 'Interpreter', 'latex');

end
