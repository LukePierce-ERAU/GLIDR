%% Full simulation of fall from 36.5km for DCR configurations
clear
close all

%% Next steps for DCR?
%   -   Develop Glide condition [Thomas]
%   -   Create a system that can do everything up to glide condition in one
%       main script [Kaylee}
%   -   impliment full 3DOF with theta-AoA relationship +++++++++++++++++++
%   -   Recreate in simulink [PDR task]
%   -   cm calc [josiah]
%   -



%% GEOMETRY ASSUMPTION SECTIONS AND GLOBAL VARIABLES ++++++++++++++++++++++


%config = 1; % Change to run sims on each config
% 1 ==== Blended Wing Body
% 2 ==== Rogallo Wing
% 3 ==== Bullet Bill

DESIGN.m = 7;
    DESIGN.c = 0.5; % Set number
    DESIGN.b = 2;
    DESIGN.S = 0.5; % Very crude guess based on prelim CAD work
    DESIGN.AR = DESIGN.b.^2./DESIGN.S; % Hard to guess right now
    DESIGN.e = 4.61.*(1-.045.*DESIGN.AR.^.68).*(cos(pi/6))-3.1; % Low Estimate after looking at sources VHANGE EQ HERE
    DESIGN.cLalpha = 180/pi.* 0.1;
    DESIGN.c = 1.2; % Set number
    DESIGN.S_ratio = 2.5; % Made from similar aircraft; original eqation: S_wet/S_ref
    DESIGN.M_arm = 0.1;
    DESIGN.V_des = -200;
    DESIGN.c_HT = 0.2;
    DESIGN.l_ht = 1/3;
    DESIGN.M_armH = .48;
    DESIGN.SH = (DESIGN.c_HT.*DESIGN.c.*DESIGN.S)/DESIGN.l_ht;
    DESIGN.g = 9.81;







% Freefall until desired q
options = odeset('Events', @y1_free);
place_free = @(t,y) funfree(t,y,DESIGN);
[t,y,te,ye,ei] = ode45(place_free, linspace(0,100,1000), [36000;0], options); % need to impliment aoa for linear region

% free(:,:,ii) = {t;y(:,1);y(:,2)};
% free_out(:,:,ii) = [te,ye,ei];

global pull_init
pull_init = te;

% Pullout @ desired q
options = odeset('Events', @y1_pull);
place_pull = @(t,y) funpull(t,y,DESIGN);
[t,y,te,ye,ei] = ode45(place_pull, [0 100], [0;ye(1,1);ye(1,2);-pi/2;0], options); % need to impliment theta and Q

% % pull(:,:,ii) = {t;y(:,1);y(:,2);y(:,3);y(:,4),y(:,5)};
% pull_out(:,:,ii) = [te,ye,ei];
% Glide @ V for L/D max


% Landing (Bleed V for safe landing)



%% Cool Plots

% Plotting tools
figure
plot(y(:,1),y(:,2))



%% Functions for simulations


% Freefall section ode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function dydt = funfree(t,y,DESIGN)


% Get atmospheric properties at current height
    [T, a, P, rho, nu, mu] = atmosisa(y(1), 'extended', true);

    % Drag Buildup for changing Re and flight condition
    Re = rho.*DESIGN.c.*y(2)./mu;

cL = 0;
[cD0,cDi] = drag(Re,cL,DESIGN);
cD = cD0 + cDi;

drag_force = 0.5 .* cD .* DESIGN.S .* rho .* y(2).^2;

    % Calculate the drag force and acceleration due to gravity
     % Drag force (assumes velocity is y(2))
    dydt = [y(2); -DESIGN.g + drag_force ./ DESIGN.m]; % [dy/dt, dv/dt]
end

function [position,isterminal,direction] = y1_free(t,y)

position = y(2) + 200; % The value that we want to be zero (altitude)
    isterminal = 1;  % Halt integration 
    direction = 0;   % The zero can be approached from either direction
end



% Pullout section ode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function dydt = funpull(t,y,DESIGN)



% Get atmospheric properties at current height
    [T, a, P, rho, nu, mu] = atmosisa(y(2), 'extended', true);

    % Drag Buildup for changing Re and flight condition
    Re = rho.*DESIGN.c.*y(3)./mu;

    % TO BE FURTHER EVALUATED++++++++++++++++++++++++++++++++++++++++
    % Estimation based on flaps to to create negative lift 
    global pull_init

    if t == pull_init
        cL = 0;
    else
        cL = 0.5;
    end
   
    
[cD0,cDi] = drag(Re,cL,DESIGN);
cD = cD0 + cDi;
    
drag_force = 0.5 .* rho .* y(3)^2 .* DESIGN.S .* DESIGN.S_ratio .* cD;
lift_force = 0.5 .* rho .* y(3)^2 .* DESIGN.S .* cL;
   
    
% Q assuming aircraft is point mass. Can be updated with inertia. Replace m
% with Izz to do so. ++++++++++++++++++++++++++++++++++++++++++++++++++++

% Assumes that aero center is on same waterline as cg.++++++++++++++++++++

Q_dot = lift_force./DESIGN.m.*(DESIGN.M_arm);
    
 % Calculate the drag force and acceleration due to gravity
     % Drag force (assumes velocity is y(2))
    dydt = [y(2).*cos(y(4)); y(2).*sin(y(4)); -DESIGN.g.*sin(y(4)) + drag_force / DESIGN.m; y(5); Q_dot]; % [dx/dt, dy/dt, dv/dt, dtheta/dt, dQ/dt]
end

function [position,isterminal,direction] = y1_pull(t,y)
    position = y(4); % The value that we want to be zero (altitude)
    isterminal = 1;  % Halt integration 
    direction = 0;   % The zero can be approached from either direction
end