function setup()
% setup - Adds the demo folders to the MATLAB path and checks for the Control System Toolbox
%
% Inputs:
%    none
%
% Outputs:
%    none

demoRoot = fileparts(mfilename('fullpath'));
for folder = ["lqrDesign", "systems", "simulation", "plotting"]
    addpath(genpath(fullfile(demoRoot, folder)));
end

if ~exist('idare', 'file')
    error('setup:noControlToolbox', ...
        'idare not found: the demo requires the Control System Toolbox.');
end

end
