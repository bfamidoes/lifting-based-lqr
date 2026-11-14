function plotFormationResults(xLifted, xStepwise, p)
% plotFormationResults - Chaser position relative to the chief for both feedback laws
%
% Position components against time, and the radial/along-track and radial/cross-track planes.
%
% Inputs:
%    xLifted, xStepwise - LVLH states (km, km/s), from simulateFormationNonlinear
%    p                  - Steps per orbit
%
% Outputs:
%    none (two figures)

orbit = (0:size(xLifted, 2) - 1) / p;
axisNames = ["Radial", "Along-track", "Cross-track"];

figure('Name', 'Tracking errors');
for i = 1:3
    subplot(3, 1, i);
    plot(orbit, xStepwise(i, :), 'b', orbit, xLifted(i, :), 'r');
    grid on;
    ylabel(axisNames(i) + " (km)", 'Interpreter', 'latex');
end
legend('Stepwise feedback', 'Lifted feedback', 'Interpreter', 'latex', 'Location', 'best');
xlabel('Orbit', 'Interpreter', 'latex');

figure('Name', 'Tracking error slices');
for i = 1:2
    subplot(1, 2, i);
    plot(xStepwise(i + 1, :), xStepwise(1, :), 'b', xLifted(i + 1, :), xLifted(1, :), 'r');
    hold on;
    plot(xStepwise(i + 1, 1), xStepwise(1, 1), 'k.', 'MarkerSize', 20);
    plot(0, 0, 'g.', 'MarkerSize', 15);
    grid on;
    xlabel(axisNames(i + 1) + " (km)", 'Interpreter', 'latex');
    ylabel('Radial (km)', 'Interpreter', 'latex');
end
legend('Stepwise feedback', 'Lifted feedback', 'Initial position', 'Chief', 'Interpreter', 'latex', 'Location', 'best');

end
