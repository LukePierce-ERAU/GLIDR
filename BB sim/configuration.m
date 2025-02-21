function [DESIGN] = configuration(config)
% Configuration selection and defining geometry for each
%   Configurations are chosen and for all stages of flight are organized
%   
%   First set of data is configuration specific. What follows if general constants true for all designs

%   Configuration specific array order is as follows:
%   1 =    m       --->   mass [kg]
%   2 =    AR      --->   Aspect ratio []
%   3 =    e       --->   Oswold's Efficiency'e Factor []
%   4 =    cLalpha --->   3D Lift curve slope [1/rad]
%   5 =    c       --->   chord length (wing) [m]
%   6 =    S       --->   Wing reference area [m^2]
%   7 =    S_ratio --->   Wing reference area to total wetted area [m^2]
%   8 =    M_arm   --->   Center of Lift moment arm from CG [m]
%   9 =    V_des   --->   Desired pull-up trim speed [m/s]

%   General constants that follow:
%   1 =    g       --->   Gravity [m/s^2]
%   2 =    eta_h   --->   Downwash factor



if config == 1
    % Geometry assumptions for Blended Wing Body

    DESIGN.m = 8;
    DESIGN.c = 1.2; % Set number
    DESIGN.b = 0.45;
    DESIGN.S = (linspace(1,2,10))'; % Very crude guess based on prelim CAD work
    DESIGN.AR = DESIGN.b.^2./DESIGN.S; % Hard to guess right now
    DESIGN.e = 0.35; % Low Estimate after looking at sources VHANGE EQ HERE
    DESIGN.cLalpha = 180/pi.* 0.01;
    DESIGN.S_ratio = 2.5; % Made from similar aircraft; original eqation: S_wet/S_ref
    DESIGN.M_arm = 0.1;
    DESIGN.V_des = -200;
    

elseif config == 2
    % Geometry assumptions for Rogallo Wing
    DESIGN.m = 10;
    DESIGN.c = 0.5; % Set number
    DESIGN.b = 2;
    DESIGN.S = (linspace(0.5,2,10))'; % Very crude guess based on prelim CAD work
    DESIGN.AR = DESIGN.b.^2./DESIGN.S; % Hard to guess right now
    DESIGN.e = 4.61.*(1-.045.*DESIGN.AR.^.68).*(cos(pi/6))-3.1; % Low Estimate after looking at sources VHANGE EQ HERE
    DESIGN.cLalpha = 180/pi.* linspace(0.1,0.07,10);
    DESIGN.c = 1.2; % Set number
    DESIGN.S_ratio = 2.5; % Made from similar aircraft; original eqation: S_wet/S_ref
    DESIGN.M_arm = 0.1;
    DESIGN.V_des = -200;
    DESIGN.c_HT = 0.2;
    DESIGN.l_ht = 1/3;
    DESIGN.M_armH = .48;
    DESIGN.SH = (DESIGN.c_HT.*DESIGN.c.*DESIGN.S)/DESIGN.l_ht;

    DESIGN.cLaH = .707;
    DESIGN.nH = 1;
    DESIGN.SH = .6;
    DESIGN.de = 2/3;
    DESIGN.M_armH = .048;
    DESIGN.CmAC = 0 ; % 
    DESIGN.CLO = .3 ; % 
    DESIGN.tauE = .5;


elseif config == 3
    % Geometry assumptions for Bullet Bill
    DESIGN.m = 5;      % Mass (kg)
    DESIGN.S = linspace(1,3,3);
    DESIGN.c = 1.2;      % Chord Length (m)
    DESIGN.b = 0.45;     % Wingspan (m)
    DESIGN.AR = DESIGN.b.^2 ./ DESIGN.S; % Aspect Ratio
    DESIGN.e = 0.4;     % Oswald Efficiency Factor
    DESIGN.cLalpha = (180/pi * 0.01);
    % DESIGN.cLalpha = 180/pi * 0.01; % Lift Curve Slope
    DESIGN.S_ratio = 2.5;
    %DESIGN.M_arm = 0.1;   % Center of Lift Moment Arm from CG (m)
    %DESIGN.V_des = -200;  % Desired Pull-Up Trim Speed (m/s)
    % DESIGN = [1, 2, 3]; % All assumptions make an array
end
     DESIGN.g = 9.81;
 % DESIGN.cLalpha = 5.7; % Lift curve slope

end

     % DESIGN.g = 9.81;

     % end