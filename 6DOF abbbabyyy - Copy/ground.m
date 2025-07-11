function [value, isterminal, direction] = ground(t, y, DESIGN)
 
    % Altitude (assuming z-down): target ground altitude
    altitude = y(3);  % If s = [x; y; z], and z increases downward

    % Event 1: Deploy paraglider at a preset altitude
    break_alt = DESIGN.break_alt;
    value = altitude + break_alt; % Variables already have opposite signs
    isterminal = 1;      % Stop integration at landing
    direction = 0;      % Either direction termination
end