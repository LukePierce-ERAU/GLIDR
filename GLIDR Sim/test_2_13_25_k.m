%% =========================================================================
%  FLOATing DRAGON - Freefall & Pullout Simulation (MATLAB Script)
%  This script simulates the descent of a vehicle from 36.5 km altitude
%  It models both:
%   - Freefall (uncontrolled descent)
%   - Pullout (controlled deceleration using lift & drag)
%  The script uses ODE45 to solve the equations of motion.
% =========================================================================
clear
close all
%% ======================== GLOBAL PARAMETERS ==============================
% Configurations:
%   1 - Blended Wing Body
%   2 - Rogallo Wing
%   3 - Bullet Bill
config = 2;  % Change this value to select a configuration
% Load design parameters from configuration file
DESIGN = configuration(config);  
end_sim = size(DESIGN.S,1);  % Number of simulation steps
%% ====================== INITIAL CONDITIONS ===============================
% State vector: [altitude (m), velocity (m/s), theta (rad), angular velocity (rad/s)]
y0 = [36000;  % Initial altitude (m)
      0;      % Initial velocity (m/s) - starts from rest
      -pi/2;  % Initial orientation angle (radians) - downward
      0];     % Initial angular velocity (rad/s)
global cL_W cL_H
cL_W = 0.1;
cL_H = 0.5;
%% ===================== ODE45 SIMULATION SETUP ============================
% The simulation runs until the vehicle reaches 100m altitude.
% This ensures we only model freefall & pullout (glide phase handled separately).
for ii = 1:1   %
options = odeset('Events', @y1_event); % Stop condition for ODE45
% Solve the equations of motion using ODE45
[t, y] = ode45(@(t, y) combinedEoM(t, y, DESIGN, ii, config), [0 5000], y0, options);
end
%% ======================= PLOTTING RESULTS ================================
% 2D Flight Path Visualization
figure
plot(y(:,3), y(:,1), 'b'); % Theta vs Altitude
xlabel('Theta (rad)');
ylabel('Altitude (m)');
title('Combined Freefall & Pullout');
%% ================== EQUATIONS OF MOTION FUNCTION =========================
function dydt = combinedEoM(t, y, DESIGN, ii, config)
    % --------------------------------------------------------------
    % Function: combinedEoM
    % This function computes the equations of motion for freefall & pullout.
    % Inputs:
    %   - t: Current time (s)
    %   - y: State vector [alt, vel, theta, omega]
    %   - DESIGN: Struct containing vehicle parameters
    %   - ii: Configuration index
    % Outputs:
    %   - dydt: Time derivatives of state variables
    % --------------------------------------------------------------
    global cL_W cL_H
    % Extract state variables
    alt = y(1);      % Altitude (m)
    vel = y(2);      % Velocity (m/s)
    theta = y(3);    % Orientation Angle (radians)
    omega = y(4);    % Angular Velocity (rad/s)
    % ================== ATMOSPHERIC MODEL ===================
    % Compute air properties based on altitude
    [T, a, P, rho, nu, mu] = atmosisa(alt, 'extended', true);
    % ================== Aero force calculations  =======================
    % Compute Reynolds number (flow properties)
    Re = rho .* DESIGN.c .* abs(vel) ./ mu;
    % Set default Lift Coefficient (cL) and adjust in pullout
        %cL_alpha = 0.1; % Sensitivity of Lift to Theta
        % cL = cL_0 + cL_alpha * theta; % Adjust cL dynamically +++++++++++
    % Compute drag and lift coefficients
    [cD0, cDi] = drag(Re, cL_W, DESIGN, ii);
    cD = (cD0 + cDi)*1.3;
    
    % [cMo, cMalpha, cMiH, cMdE] = Cm(DESIGN,cL);
    % 
    % cM = 
    % Compute Aero Forces
    drag_force = 0.5 * cD * DESIGN.S(ii) * rho * vel^2;
    lift_force_W = 0.5 * cL_W * DESIGN.S(ii) * rho * vel^2;
    lift_force_H = 0.5 * cL_H * DESIGN.SH(ii) * rho * vel^2;
    % aero_moment = 0.5 * cM * DESIGN.S(ii) * rho * vel^2;
   

if theta > -pi/4
    
% % ================== COMPUTE CALCULATED VELOCITY ===================
%         % This is the expected velocity based on aerodynamic forces
%         v_calc = sqrt(2 * DESIGN.m * DESIGN.g / (rho * DESIGN.S(ii)) ...
%             * sqrt((1 / (pi * DESIGN.e * DESIGN.AR(ii))) / (3 * cD0)));
%         % ================== ADJUST DRAG IF SIMULATED VELOCITY > CALCULATED ===================
%         % If the simulated velocity is higher than the expected velocity, increase drag
% 
%         if v_calc > vel
% 
%             cL_W = cL_W *1.20;
%             cL_H = cL_H *1.20;
%         else
%             cL_W = cL_W /1.20;
%             cL_H = cL_H /1.20;
%         end

        moment = lift_force_W*DESIGN.M_arm-lift_force_H*DESIGN.M_armH;

        if theta > -pi/12 && theta < pi/12
            if moment > 0
                cL_W = cL_W /1.20;
                cL_H = cL_H *1.20;
            else
                cL_W = cL_W *1.20;
                cL_H = cL_H /1.20;
            end
        elseif theta > 0 
            cL_W = cL_W *1.20;
                cL_H = cL_H /1.20;
        end



    % BAD LOGIC +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % if theta > -pi/4
    %     %
    %     % ================== COMPUTE CALCULATED VELOCITY ===================
    %     % This is the expected velocity based on aerodynamic forces
    %     v_calc = sqrt(2 * DESIGN.m * DESIGN.g / (rho * DESIGN.S(ii)) ...
    %         * sqrt((1 / (pi * DESIGN.e * DESIGN.AR(ii))) / (3 * cD0)));
    %     % ================== ADJUST DRAG IF SIMULATED VELOCITY > CALCULATED ===================
    %     % If the simulated velocity is higher than the expected velocity, increase drag
    % 
    %     if abs(vel) > v_calc
    %         cL_H = cL_H * 1.01;  % Increase Lift (and indirectly Drag) by 5%
    %     else 
    %         cL_H = cL_H / 1.01;
    %     end
    % 
    %     moment = lift_force_W*DESIGN.M_arm-lift_force_H*DESIGN.M_armH
    % 
    %     if theta > 0 && abs(moment) >  2
    % 
    %         if moment > 0
    %         cL_W = cL_W / 1.10 % Increase Lift (and indirectly Drag) by 5%
    %         else
    %         cL_W = cL_W * 1.10
    %         end
    % 
    % 
    %     end
    % 
    % 
    % 
    %     %
    % end
end

    % ================== COMPUTE EQUATIONS OF MOTION ===================
    % Vertical Acceleration (gravity + drag)
    accel = -DESIGN.g + drag_force / DESIGN.m;
    
    % Angular Acceleration (depends on lift)
    alpha = (lift_force_W-lift_force_H) / DESIGN.m * DESIGN.M_arm;
    % ================== RETURN TIME DERIVATIVES ===================
    dydt = [vel;      % dy/dt (Altitude change)
            accel;    % dv/dt (Velocity change)
            omega;    % dtheta/dt (Orientation change)
            alpha];   % domega/dt (Angular velocity change)
end
%% ====================== EVENT FUNCTION (STOP SIMULATION) =========================
function [position,isterminal,direction] = y1_event(t,y)
    % --------------------------------------------------------------
    % Function: y1_event
    % Stops the simulation when altitude reaches 100m.
    % --------------------------------------------------------------
    position = y(1) - 100; % Stop when reaching 100m altitude
    isterminal = 1;  % Halt integration
    direction = 0;   % The zero can be approached from either direction
end