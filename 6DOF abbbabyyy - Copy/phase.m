function [value, isterminal, direction] = phase(t, y, DESIGN)

    % Altitude (assuming z-down): target ground altitude
    altitude = y(3);  % If s = [x; y; z], and z increases downward

    % Event 1: Deploy paraglider at a preset altitude
    glidr_alt = DESIGN.glidr_alt;
    value = altitude + glidr_alt; % Variables already have opposite signs
    isterminal = 1;      % Stop integration to switch
    direction = 0;      % Either direction termination

end