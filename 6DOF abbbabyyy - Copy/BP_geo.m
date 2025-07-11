function [DESIGN] = BP_geo()
% Design definition ~~~~~~~~~~~~~~~~ 
% ALL NUMBERS TO BE CHECKED AND UPDATED WEEK OF 8/25/25!!!
% Body cordinate system is defined as x out of the bottom of the crumple
% zone, y out of the right side of the paraglider when "behind" trailing
% edge. Z faces opposite direction of the front airbrake when deployed.
% (Defined for geometry definitions)


DESIGN.m = 5; % kg
DESIGN.FF_I11 = 03.3;
DESIGN.FF_I13 = 00.03;
DESIGN.FF_I22  = 03.4;
DESIGN.FF_I33 = 01.90;
DESIGN.glidr_alt = 18000;
DESIGN.break_alt = 1700;

% Airbrake geometry and information
DESIGN.ab_S = 1; % [m^2] Reference area of one airbrake
DESIGN.ab_count = 3; % Number of airbrakes in design
DESIGN.ab_separation = 120; % [degree] Airbrake separation angle between airbrakes;
DESIGN.ab_arcdeg = 120; % [degree] Airbrake span angle in degrees
DESIGN.ab_b = 0.9; % [m]
DESIGN.ab_c = 0.15; % [m]
DESIGN.ab_t = 0.07; % [m] thickness of airbrake
DESIGN.ab_r = 0.15; % [m] Radius of the airbrake curve
DESIGN.ab_arc = DESIGN.ab_arcdeg*pi/180 * DESIGN.ab_r; % [m] arc length of the airbrake
DESIGN.ab_Swet = DESIGN.ab_arc*2*DESIGN.ab_b; % [m^2] wetted area for a single airbrake
DESIGN.ab_def_max = 10*pi/180; % [rad] max airbrake deflection in radians

DESIGN.ab_cLalpha = 2*pi; % [rad^-1] This is a pure guess


% Geometry of fuselage
DESIGN.f_h = 0.9; % [m] Height of fuelage
DESIGN.f_d = 0.35; % [m] Diameter of fuselage
DESIGN.f_Sref = pi/4 * DESIGN.f_d^2;
DESIGN.f_Swet = 2*DESIGN.f_Sref + DESIGN.f_h * pi*DESIGN.f_d; % [m^2] Wetted area of fuselage

% Geometry of airbrake mounting system
DESIGN.ms_l = 0.05; % [m] Length of protruding braket for airbrake mounting

% Geometry of glidr
DESIGN.glidr_S_ref = 1.57; % [m^2] Projected Reference area of glidr
DESIGN.glidr_S_arc = 1.89; % [m^2] Flat reference area of glidr
DESIGN.glidr_S_wet = 4.25 * DESIGN.glidr_S_arc; % [m^2] Wetted area of glidr
% DESIGN.glidr_b = ;
% DESIGN.glidr_c = 

% Full system body geometry
DESIGN.FF_sb_b = 2*(DESIGN.f_d/2 + DESIGN.ms_l + DESIGN.ab_b)*cosd(30);                                            % [m] Span of body when fuselage is deployed. NOTE: 30 degrees is dependent on number of airbrakes
DESIGN.FF_sb_c = (DESIGN.f_d/2 + DESIGN.ms_l + DESIGN.ab_b) + (DESIGN.f_d/2 + DESIGN.ms_l + DESIGN.ab_b)*sind(30); % [m] chord along x body axis with airbrakes deployed
DESIGN.FF_sb_Sref = 3*DESIGN.ab_S + pi/4*DESIGN.f_d^2; % [m^2] total reference area of body in freefall

% Target location (PLEASE CHANGE THIS SECTION BRUCEEE) only angle part tho
DESIGN.target_psi = 20  *pi/180; % To be based on chords later...


end