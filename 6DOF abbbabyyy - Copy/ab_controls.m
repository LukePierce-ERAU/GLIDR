function [delta_ab] = ab_controls(DESIGN,psi)
%UNTITLED Summary of this function goes here
%   Proportional Controls for Airbrake controls

K_ab = 1.5; % Gain value

error = DESIGN.target_psi - psi;

delta_ab = K_ab * error;

if delta_ab > DESIGN.ab_def_max
    delta_ab = DESIGN.ab_def_max;
else

end