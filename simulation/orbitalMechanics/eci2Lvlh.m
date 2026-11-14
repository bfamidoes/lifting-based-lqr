function [pRel, vRel, eciToLvlhDcm] = eci2Lvlh(rChief, vChief, rChaser, vChaser)
% eci2Lvlh - Chaser state relative to the chief, in the chief's LVLH frame
%
% LVLH: x radial, y along-track, z orbit normal; the frame rotates at |h|/|r|^2 under two-body motion.
%
% Inputs:
%    rChief, vChief   - Chief position and velocity in ECI (m, m/s)
%    rChaser, vChaser - Chaser position and velocity in ECI (m, m/s)
%
% Outputs:
%    pRel, vRel   - Relative position and velocity in LVLH (m, m/s)
%    eciToLvlhDcm - Rotation matrix from ECI to LVLH

h = cross(rChief, vChief);
xHat = rChief / norm(rChief);
zHat = h / norm(h);
yHat = cross(zHat, xHat);
eciToLvlhDcm = [xHat'; yHat'; zHat'];

omega = [0; 0; norm(h) / norm(rChief)^2];
pRel = eciToLvlhDcm * (rChaser - rChief);
vRel = eciToLvlhDcm * (vChaser - vChief) - cross(omega, pRel);

end
