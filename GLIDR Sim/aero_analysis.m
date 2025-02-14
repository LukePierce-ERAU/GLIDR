config = 1;


% Aero analysis Assumptions for BL vehicle
    % m = 8;
    % b = 0.45;
    % S = linspace(1,2,10); % Very crude guess based on prelim CAD work
    % AR = b^2./S;
    % S_ratio = 2.5; % Made from similar aircraft; original eqation: S_wet/S_ref
    % e = 0.35;
    % cLalpha = 1./(pi.*e.*AR);
    

% Paraglider (fullsize)
    % m = 6;
    % AR = [4.5 4.6 4.7 4.8 4.9 5 5.1 5.2 5.3 5.4 5.5];
    % e = 1.78.*(1-.045.*AR.^.68)-.64; % Oswald Efficiency
    % S = linspace(21,24.5); % m^2, where did you get this from
    % cLalpha = 1/(pi*e*AR);

% Rogallo Wing (fullsize)
    m = 7;
    b = 2;
    S = linspace(0.5,1,10); % m^2
    AR = b.^2./S;
    e = 4.61.*(1-.045.*AR.^.68).*(cos(pi/6))-3.1; % Oswald Efficiency
    cLalpha = 1./(pi.*e.*AR);