%% Full simulation of fall from 36.5km for DCR configurations
clear
close all
 
%% Next steps for DCR?
%   -   Combine Freefall & Pullout into a single ODE45 simulation ✅ (Done)
%   -   Ensure outputs can handle theta/alpha relationship when ready ✅
%   -   Remove separate Freefall & Pullout ODE45 calls ✅
%   -   Leave Glide Phase for later integration ✅ (Glide Removed)
 
%% GEOMETRY ASSUMPTION SECTIONS AND GLOBAL VARIABLES ++++++++++++++++++++++
 
config = 1; % Change to run different configurations
% 1 ==== Blended Wing Body
% 2 ==== Rogallo Wing
% 3 ==== Bullet Bill
 
DESIGN = configuration(config);  % Load specific configuration parameters
end_sim = size(DESIGN.S,1);  % Get total simulation steps
 
%% INITIAL CONDITIONS
y0 = [36000;  % Initial altitude (m)
      0;      % Initial velocity (m/s)
      -pi/2;      % Initial theta (radians)
      0];     % Initial angular velocity (rad/s)
 
%% MAIN SIMULATION LOOP
for ii = 1:1:end_sim
 
    % 🟢 FREEFALL & PULLOUT PHASE: Single ODE45 simulation
    options = odeset('Events', @y1_event); % Event-based stopping condition
    place_combined = @(t,y) combinedEoM(t, y, DESIGN, ii);
 
    [t, y] = ode45(place_combined, [0 500], y0, options);
 
    % Store Simulation Results
    sim_data(:,:,ii) = {t, y};
 
end
 
%% Cool Plots
% 🟣 2D Flight Path Visualization
figure
plot(y(:,3), y(:,1), 'b'); 
xlabel('Theta'); 
ylabel('Altitude (m)');
title('Combined Freefall & Pullout');
 
%% FUNCTIONS FOR SIMULATION
 
% 🟢 COMBINED EQUATIONS OF MOTION FOR FREEFALL + PULLOUT ~~~~~~~~~~~~~~~~
function dydt = combinedEoM(t, y, DESIGN, ii)
    % Extract state variables
    alt = y(1);      % Altitude
    vel = y(2);      % Velocity
    theta = y(3);    % Orientation Angle (theta)
    omega = y(4);    % Angular Velocity
 
    % Get atmospheric properties at current height
    [T, a, P, rho, nu, mu] = atmosisa(alt, 'extended', true);
 
    % Compute Drag
    
   Re = rho .* DESIGN.c .* abs(vel) ./ mu;
    
   % Compute Lift (Activated During Pullout)
    
   
   
   if theta > -pi/4  % If pullout is triggered
        cL =   % Placeholder Lift Coefficient for BL
    else
        cL = 0.5
    end

    [cD0,cDi] = drag(Re, cL, DESIGN, ii);
    cD = cD0 + cDi;
    drag_force = 0.5 * cD * DESIGN.S(ii) * rho * vel^2;
  
    lift_force = 0.5 * cL * DESIGN.S(ii) * rho * vel^2;
 
    % Compute Acceleration
    accel = -DESIGN.g + drag_force / DESIGN.m;  % Vertical Acceleration
    alpha = lift_force / DESIGN.m * DESIGN.M_arm;  % Angular Acceleration
 
    % Return dydt
    dydt = [vel; accel; omega; alpha]; % [dy/dt, dv/dt, dtheta/dt, domega/dt]
end
 
% 🔵 EVENT FUNCTION TO STOP SIMULATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [position,isterminal,direction] = y1_event(t,y)
    position = y(1) - 100; % Stop when reaching 100m altitude
    isterminal = 1;  % Halt integration 
    direction = 0;   % The zero can be approached from either direction
end