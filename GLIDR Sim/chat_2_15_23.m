

    % % Aircraft parameters
    % m = 5000;                % Mass of the aircraft (kg)
    % S = 30;                  % Wing reference area (m²)
    % C_L = 0.5;               % Coefficient of lift (assumed for glide)
    % C_D = 0.02;              % Coefficient of drag (assumed for glide)
    % rho = 1.225;             % Air density at sea level (kg/m³)
    % g = 9.81;                % Gravitational acceleration (m/s²)
    % 
    % % Flight conditions
    % altitude = 1000;         % Altitude (m) (used for adjusting air density if needed)
    % 
    % % Compute weight
    % W = m * g;               % Weight of the aircraft (N)
    % 
    % % Calculate steady glide velocity (assuming a balance between lift and weight)
    % % For glide, lift equals weight, and the descent path angle is determined by D/L
    % % Find the glide speed based on the balance of drag and lift
    % % The glide ratio is approximately: L/D = V^2 / (m * g)
    % glide_ratio = (2 * W) / (rho * S * (C_L / C_D));  % Glide ratio
    % 
    % % We can find the steady-state glide velocity (V) using the glide ratio
    % V_glide = sqrt((2 * W) / (rho * S * (C_L / C_D))); % Glide velocity (m/s)
    % 
    % % Compute lift and drag at the glide velocity
    % L = 0.5 * rho * V_glide^2 * C_L * S;  % Lift force (N)
    % D = 0.5 * rho * V_glide^2 * C_D * S;  % Drag force (N)
    % 
    % % Calculate glide path angle (gamma)
    % gamma = atan(D / L);  % Flight path angle (radians)
    % 
    % % Output results
    % fprintf('Steady-State Glide Simulation Results:\n');
    % fprintf('Glide Speed (V): %.2f m/s\n', V_glide);
    % fprintf('Lift (L): %.2f N\n', L);
    % fprintf('Drag (D): %.2f N\n', D);
    % fprintf('Glide Path Angle (gamma): %.2f degrees\n', rad2deg(gamma));
    % 
    % % Plotting the glide path and forces
    % figure;
    % subplot(1,2,1);
    % plot([0, 1000], [0, -1000*tan(gamma)], 'LineWidth', 2);
    % title('Steady-State Glide Path');
    % xlabel('Distance (m)');
    % ylabel('Altitude (m)');
    % grid on;
    % 
    % subplot(1,2,2);
    % bar([L, D, W]);
    % set(gca, 'xticklabel', {'Lift', 'Drag', 'Weight'});
    % title('Forces in Steady Glide');
    % ylabel('Force (N)');


    % % Aircraft parameters
    % m = 500;                % Mass of the aircraft (kg)
    % S = 30;                  % Wing reference area (m²)
    % C_L0 = 0.3;              % Lift coefficient at zero angle of attack (dimensionless)
    % C_L_alpha = 0.1;         % Lift curve slope (per radian) (dimensionless)
    % C_D0 = 0.02;             % Parasitic drag coefficient (dimensionless)
    % k = 0.04;                % Induced drag factor (dimensionless)
    % rho = 1.225;             % Air density at sea level (kg/m³)
    % g = 9.81;                % Gravitational acceleration (m/s²)
    % 
    % % Flight conditions
    % altitude = 1000;         % Altitude (m) (used for adjusting air density if needed)
    % 
    % % Compute weight
    % W = m * g;               % Weight of the aircraft (N)
    % 
    % % Glide velocity initial guess
    % V_glide = 50;            % Initial guess for glide velocity (m/s)
    % tolerance = 0.01;        % Convergence tolerance for velocity
    % 
    % % Iterate to find the steady-state glide velocity where lift equals weight
    % iteration = 0;
    % V_new = V_glide;
    % while true
    %     % Compute angle of attack (assuming glide angle is small, use L/D balance)
    %     alpha = atan(W / (0.5 * rho * V_new^2 * S * C_L0));  % Approximate angle of attack
    % 
    %     % Compute changing lift and drag coefficients based on angle of attack
    %     C_L = C_L0 + C_L_alpha * alpha;
    %     C_D = C_D0 + k * C_L^2;
    % 
    %     % Calculate steady glide conditions
    %     L = 0.5 * rho * V_new^2 * C_L * S;  % Lift force (N)
    %     D = 0.5 * rho * V_new^2 * C_D * S;  % Drag force (N)
    % 
    %     % Balance forces for steady glide (Lift = Weight)
    %     velocity_required = sqrt((2 * W) / (rho * S * C_L));
    % 
    %     % Check for convergence (if velocity changes are small enough)
    %     if abs(V_new - velocity_required) < tolerance
    %         break;
    %     end
    % 
    %     % Update glide velocity for the next iteration
    %     V_new = velocity_required;
    %     iteration = iteration + 1;
    %     if iteration > 100
    %         disp('Warning: Convergence may not be achieved');
    %         break;
    %     end
    % end
    % 
    % % Output results
    % fprintf('Steady-State Glide Simulation Results with Variable Coefficients:\n');
    % fprintf('Glide Speed (V): %.2f m/s\n', V_new);
    % fprintf('Lift (L): %.2f N\n', L);
    % fprintf('Drag (D): %.2f N\n', D);
    % fprintf('Angle of Attack (alpha): %.2f degrees\n', rad2deg(alpha));
    % 
    % % Plotting the glide path and forces
    % figure;
    % subplot(1,2,1);
    % plot([0, 1000], [0, -1000*tan(alpha)], 'LineWidth', 2);
    % title('Steady-State Glide Path');
    % xlabel('Distance (m)');
    % ylabel('Altitude (m)');
    % grid on;
    % 
    % subplot(1,2,2);
    % bar([L, D, W]);
    % set(gca, 'xticklabel', {'Lift', 'Drag', 'Weight'});
    % title('Forces in Steady Glide with Varying Coefficients');
    % ylabel('Force (N)');



    % % Aircraft parameters
    % m = 5000;                % Mass of the aircraft (kg)
    % S = 30;                  % Wing reference area (m²)
    % C_L0 = 0.3;              % Lift coefficient at zero angle of attack (dimensionless)
    % C_L_alpha = 0.1;         % Lift curve slope (per radian) (dimensionless)
    % C_D0 = 0.02;             % Parasitic drag coefficient (dimensionless)
    % k = 0.04;                % Induced drag factor (dimensionless)
    % rho = 1.225;             % Air density at sea level (kg/m³)
    % g = 9.81;                % Gravitational acceleration (m/s²)
    % 
    % % Flight conditions
    % altitude = 1000;         % Altitude (m) (used for adjusting air density if needed)
    % 
    % % Compute weight
    % W = m * g;               % Weight of the aircraft (N)
    % 
    % % Glide velocity initial guess
    % V_glide = 50;            % Initial guess for glide velocity (m/s)
    % tolerance = 0.01;        % Convergence tolerance for velocity
    % 
    % % Iterate to find the steady-state glide velocity where lift equals weight
    % iteration = 0;
    % V_new = V_glide;
    % while true
    %     % Compute changing lift and drag coefficients based on the glide velocity
    %     % Angle of attack approximation from glide velocity
    %     alpha = atan(W / (0.5 * rho * V_new^2 * S * C_L0));  % Approximate angle of attack
    % 
    %     % Compute changing lift and drag coefficients based on angle of attack
    %     C_L = C_L0 + C_L_alpha * alpha;
    %     C_D = C_D0 + k * C_L^2;
    % 
    %     % Calculate steady glide conditions
    %     L = 0.5 * rho * V_new^2 * C_L * S;  % Lift force (N)
    %     D = 0.5 * rho * V_new^2 * C_D * S;  % Drag force (N)
    % 
    %     % Balance forces for steady glide (Lift = Weight)
    %     velocity_required = sqrt((2 * W) / (rho * S * C_L));  % Required velocity for steady glide
    % 
    %     % Check for convergence (if velocity changes are small enough)
    %     if abs(V_new - velocity_required) < tolerance
    %         break;
    %     end
    % 
    %     % Update glide velocity for the next iteration
    %     V_new = velocity_required;
    %     iteration = iteration + 1;
    %     if iteration > 100
    %         disp('Warning: Convergence may not be achieved');
    %         break;
    %     end
    % end
    % 
    % % Compute the glide path angle
    % gamma = atan(D / L);  % Glide path angle (radians)
    % 
    % % Output results
    % fprintf('Steady-State Glide Simulation Results with Variable Coefficients:\n');
    % fprintf('Glide Speed (V): %.2f m/s\n', V_new);
    % fprintf('Lift (L): %.2f N\n', L);
    % fprintf('Drag (D): %.2f N\n', D);
    % fprintf('Angle of Attack (alpha): %.2f degrees\n', rad2deg(alpha));
    % fprintf('Glide Path Angle (gamma): %.2f degrees\n', rad2deg(gamma));
    % 
    % % Plotting the glide path and forces
    % figure;
    % subplot(1,2,1);
    % plot([0, 1000], [0, -1000*tan(gamma)], 'LineWidth', 2);
    % title('Steady-State Glide Path');
    % xlabel('Distance (m)');
    % ylabel('Altitude (m)');
    % grid on;
    % 
    % subplot(1,2,2);
    % bar([L, D, W]);
    % set(gca, 'xticklabel', {'Lift', 'Drag', 'Weight'});
    % title('Forces in Steady Glide with Varying Coefficients');
    % ylabel('Force (N)');
    
    
    % % Aircraft parameters
    % m = 5000;                % Mass of the aircraft (kg)
    % S = 30;                  % Wing reference area (m²)
    % C_L0 = 0.3;              % Lift coefficient at zero angle of attack (dimensionless)
    % C_L_alpha = 0.1;         % Lift curve slope (per radian) (dimensionless)
    % C_D0 = 0.02;             % Parasitic drag coefficient (dimensionless)
    % k = 0.04;                % Induced drag factor (dimensionless)
    % rho = 1.225;             % Air density at sea level (kg/m³)
    % g = 9.81;                % Gravitational acceleration (m/s²)
    % 
    % % Flight conditions
    % altitude = 1000;         % Altitude (m) (used for adjusting air density if needed)
    % 
    % % Compute weight
    % W = m * g;               % Weight of the aircraft (N)
    % 
    % % Glide velocity initial guess
    % V_glide = 50;            % Initial guess for glide velocity (m/s)
    % tolerance = 0.01;        % Convergence tolerance for velocity
    % 
    % % Iterate to find the steady-state glide velocity where lift equals weight
    % iteration = 0;
    % V_new = V_glide;
    % while true
    %     % Compute changing lift and drag coefficients based on the glide velocity
    %     % Angle of attack approximation from glide velocity
    %     alpha = atan(W / (0.5 * rho * V_new^2 * S * C_L0));  % Approximate angle of attack
    % 
    %     % Compute changing lift and drag coefficients based on angle of attack
    %     C_L = C_L0 + C_L_alpha * alpha;
    %     C_D = C_D0 + k * C_L^2;
    % 
    %     % Calculate steady glide conditions
    %     L = 0.5 * rho * V_new^2 * C_L * S;  % Lift force (N)
    %     D = 0.5 * rho * V_new^2 * C_D * S;  % Drag force (N)
    % 
    %     % Balance forces for steady glide (Lift = Weight)
    %     velocity_required = sqrt((2 * W) / (rho * S * C_L));  % Required velocity for steady glide
    % 
    %     % Check for convergence (if velocity changes are small enough)
    %     if abs(V_new - velocity_required) < tolerance
    %         break;
    %     end
    % 
    %     % Update glide velocity for the next iteration
    %     V_new = velocity_required;
    %     iteration = iteration + 1;
    %     if iteration > 100
    %         disp('Warning: Convergence may not be achieved');
    %         break;
    %     end
    % end
    % 
    % % Compute the glide path angle
    % gamma = atan(D / L);  % Glide path angle (radians)
    % 
    % % Output results
    % fprintf('Steady-State Glide Simulation Results with Variable Coefficients:\n');
    % fprintf('Glide Speed (V): %.2f m/s\n', V_new);
    % fprintf('Lift (L): %.2f N\n', L);
    % fprintf('Drag (D): %.2f N\n', D);
    % fprintf('Angle of Attack (alpha): %.2f degrees\n', rad2deg(alpha));
    % fprintf('Glide Path Angle (gamma): %.2f degrees\n', rad2deg(gamma));
    % 
    % % Plotting the glide path and forces
    % figure;
    % subplot(1,2,1);
    % plot([0, 1000], [0, -1000*tan(gamma)], 'LineWidth', 2);
    % title('Steady-State Glide Path');
    % xlabel('Distance (m)');
    % ylabel('Altitude (m)');
    % grid on;
    % 
    % subplot(1,2,2);
    % bar([L, D, W]);
    % set(gca, 'xticklabel', {'Lift', 'Drag', 'Weight'});
    % title('Forces in Steady Glide with Varying Coefficients');
    % ylabel('Force (N)');

    
%% This one has better size parameters++++++++++++++++++++++++++++++++++++
    % % Aircraft parameters
    % m = 8;                % Mass of the aircraft (kg)
    % S = 2;                  % Wing reference area (m²)
    % C_L0 = 0.1;              % Lift coefficient at zero angle of attack (dimensionless)
    % C_L_alpha = 0.05;         % Lift curve slope (per radian) (dimensionless)
    % C_D0 = 0.02;             % Parasitic drag coefficient (dimensionless)
    % k = 0.12;                % Induced drag factor (dimensionless)
    % rho = 1.225;             % Air density at sea level (kg/m³)
    % g = 9.81;                % Gravitational acceleration (m/s²)
    % 
    % % Flight conditions
    % altitude = 1000;         % Altitude (m) (used for adjusting air density if needed)
    % 
    % % Compute weight
    % W = m * g;               % Weight of the aircraft (N)
    % 
    % % Glide velocity initial guess
    % V_glide = 50;            % Initial guess for glide velocity (m/s)
    % tolerance = 0.01;        % Convergence tolerance for velocity
    % 
    % % Iterate to find the steady-state glide velocity where lift equals weight
    % iteration = 0;
    % V_new = V_glide;
    % while true
    %     % Compute changing lift and drag coefficients based on the glide velocity
    %     % Angle of attack approximation from glide velocity
    %     alpha = atan(W / (0.5 * rho * V_new^2 * S * C_L0));  % Approximate angle of attack
    % 
    %     % Compute changing lift and drag coefficients based on angle of attack
    %     C_L = C_L0 + C_L_alpha * alpha;
    %     C_D = C_D0 + k * C_L^2;
    % 
    %     % Calculate steady glide conditions
    %     L = 0.5 * rho * V_new^2 * C_L * S;  % Lift force (N)
    %     D = 0.5 * rho * V_new^2 * C_D * S;  % Drag force (N)
    % 
    %     % Balance forces for steady glide (Lift = Weight)
    %     velocity_required = sqrt((2 * W) / (rho * S * C_L));  % Required velocity for steady glide
    % 
    %     % Check for convergence (if velocity changes are small enough)
    %     if abs(V_new - velocity_required) < tolerance
    %         break;
    %     end
    % 
    %     % Update glide velocity for the next iteration
    %     V_new = velocity_required;
    %     iteration = iteration + 1;
    %     if iteration > 100
    %         disp('Warning: Convergence may not be achieved');
    %         break;
    %     end
    % end
    % 
    % % Compute the glide path angle
    % gamma = atan(D / L);  % Glide path angle (radians)
    % 
    % % Angle of attack is approximately equal to the glide path angle for small angles
    % alpha = gamma;
    % 
    % % Output results
    % fprintf('Steady-State Glide Simulation Results with Variable Coefficients:\n');
    % fprintf('Glide Speed (V): %.2f m/s\n', V_new);
    % fprintf('Lift (L): %.2f N\n', L);
    % fprintf('Drag (D): %.2f N\n', D);
    % fprintf('Angle of Attack (alpha): %.2f degrees\n', rad2deg(alpha));
    % fprintf('Glide Path Angle (gamma): %.2f degrees\n', rad2deg(gamma));
    % 
    % % Plotting the glide path and forces
    % figure;
    % subplot(1,2,1);
    % plot([0, 1000], [0, -1000*tan(gamma)], 'LineWidth', 2);
    % title('Steady-State Glide Path');
    % xlabel('Distance (m)');
    % ylabel('Altitude (m)');
    % grid on;
    % 
    % subplot(1,2,2);
    % bar([L, D, W]);
    % set(gca, 'xticklabel', {'Lift', 'Drag', 'Weight'});
    % title('Forces in Steady Glide with Varying Coefficients');
    % ylabel('Force (N)');



% BAD AOA calc
    % % Aircraft parameters
    % m = 5000;                % Mass of the aircraft (kg)
    % S = 30;                  % Wing reference area (m²)
    % C_L0 = 0.3;              % Lift coefficient at zero angle of attack (dimensionless)
    % C_L_alpha = 0.1;         % Lift curve slope (per radian) (dimensionless)
    % C_D0 = 0.02;             % Parasitic drag coefficient (dimensionless)
    % k = 0.04;                % Induced drag factor (dimensionless)
    % rho0 = 1.225;            % Sea level air density (kg/m³)
    % g = 9.81;                % Gravitational acceleration (m/s²)
    % 
    % % Flight conditions (inputs)
    % altitude = 3000;         % Altitude (m)
    % temperature = 288.15 - 0.0065 * altitude; % Temperature from ISA model (Kelvin)
    % 
    % % Calculate air density at given altitude using ISA model
    % rho = rho0 * (1 - 0.0065 * altitude / 288.15) ^ 4.256;
    % 
    % % Compute weight
    % W = m * g;               % Weight of the aircraft (N)
    % 
    % % Stall angle (radians) (stall angle assumed for typical subsonic aircraft)
    % stall_angle = deg2rad(15);  
    % 
    % % Target glide velocity (optimal velocity determined earlier)
    % % Assuming the optimized value was calculated previously
    % target_velocity = sqrt((2 * W) / (rho * S * C_L));     % Example target velocity in m/s (can be calculated from previous results)
    % 
    % % Initial condition: Aircraft starts with a higher glide velocity
    % V_current = 100;           % Higher initial glide velocity (m/s)
    % tolerance = 0.01;         % Convergence tolerance for velocity
    % distance_traveled = 0;    % Distance traveled (m)
    % 
    % % Define lift coefficient limits (realistic values)
    % CL_max = 2.5;             % Maximum CL (typically around this for subsonic aircraft with flaps)
    % CL_min = 0.5;             % Minimum CL (reasonable lower bound)
    % 
    % % Iterate until the velocity is trimmed to the target
    % iteration = 0;
    % while abs(V_current - target_velocity) > tolerance
    %     % Calculate angle of attack (alpha) to match current velocity with weight
    %     alpha = atan(W / (0.5 * rho * V_current^2 * S * C_L0)); 
    % 
    %     % Calculate dynamic CL (lift coefficient) based on AoA and flap settings
    %     % Adjust CL to keep it within reasonable bounds
    %     C_L = C_L0 + C_L_alpha * alpha;
    % 
    %     % Ensure CL stays within reasonable limits (0.5 < CL < 2.5)
    %     if C_L > CL_max
    %         C_L = CL_max;
    %     elseif C_L < CL_min
    %         C_L = CL_min;
    %     end
    % 
    %     % Induced drag model (quadratic relation to lift coefficient)
    %     C_D = C_D0 + k * C_L^2;
    % 
    %     % Calculate steady glide conditions
    %     L = 0.5 * rho * V_current^2 * C_L * S;  % Lift force (N)
    %     D = 0.5 * rho * V_current^2 * C_D * S;  % Drag force (N)
    % 
    %     % Calculate distance traveled during this step
    %     distance_traveled = distance_traveled + V_current * 1;  % Assume each iteration is 1 second
    % 
    %     % Balance forces for steady glide (Lift = Weight)
    %     % Adjust AoA to reduce velocity if it's too high
    %     V_new = V_current - D / m * 1;  % Apply drag force to reduce speed (assume 1 second steps)
    % 
    %     % Update current velocity and ensure it stays above the target
    %     V_current = max(V_new, target_velocity);
    % 
    %     % Prevent angle of attack from reaching stall
    %     if alpha > stall_angle
    %         alpha = stall_angle;
    %     end
    % 
    %     % Increment iteration count
    %     iteration = iteration + 1;
    % 
    %     % Display progress
    %     if mod(iteration, 10) == 0
    %         fprintf('Iteration %d: Velocity = %.2f m/s, Distance Traveled = %.2f m, CL = %.2f, Alpha = %.2f deg\n', ...
    %             iteration, V_current, distance_traveled, C_L, rad2deg(alpha));
    %     end
    % end
    % 
    % % Final Output
    % fprintf('Final Trimming Results:\n');
    % fprintf('Final Glide Velocity: %.2f m/s\n', V_current);
    % fprintf('Total Distance Traveled: %.2f m\n', distance_traveled);
    % fprintf('Final Lift Coefficient (CL): %.2f\n', C_L);
    % fprintf('Final Angle of Attack (Alpha): %.2f degrees\n', rad2deg(alpha));
    % 
    % % Plotting the velocity reduction over time
    % figure;
    % plot(1:iteration, linspace(100, V_current, iteration), 'LineWidth', 2);
    % xlabel('Iteration');
    % ylabel('Velocity (m/s)');
    % title('Reduction in Glide Velocity Over Time');
    % grid on;





    % % Aircraft parameters
    % m = 5000;                % Mass of the aircraft (kg)
    % S = 30;                  % Wing reference area (m²)
    % C_L0 = 0.3;              % Lift coefficient at zero angle of attack (dimensionless)
    % C_L_alpha = 0.1;         % Lift curve slope (per radian) (dimensionless)
    % C_D0 = 0.02;             % Parasitic drag coefficient (dimensionless)
    % k = 0.04;                % Induced drag factor (dimensionless)
    % rho0 = 1.225;            % Sea level air density (kg/m³)
    % g = 9.81;                % Gravitational acceleration (m/s²)
    % 
    % % Flight conditions (inputs)
    % altitude = 3000;         % Altitude (m)
    % temperature = 288.15 - 0.0065 * altitude; % Temperature from ISA model (Kelvin)
    % 
    % % Calculate air density at given altitude using ISA model
    % rho = rho0 * (1 - 0.0065 * altitude / 288.15) ^ 4.256;
    % 
    % % Compute weight
    % W = m * g;               % Weight of the aircraft (N)
    % 
    % % Stall angle (radians) (stall angle assumed for typical subsonic aircraft)
    % stall_angle = deg2rad(15);  
    % 
    % % Target glide velocity (optimal velocity determined earlier)
    % % Assuming the optimized value was calculated previously
    % target_velocity = sqrt((2 * W) / (rho * S * C_L));     % Example target velocity in m/s (can be calculated from previous results)
    % 
    % % Initial condition: Aircraft starts with a higher glide velocity
    % V_current = 100;           % Higher initial glide velocity (m/s)
    % tolerance = 0.01;         % Convergence tolerance for velocity
    % distance_traveled = 0;    % Distance traveled (m)
    % 
    % % Define lift coefficient limits (realistic values)
    % CL_max = 2.5;             % Maximum CL (typically around this for subsonic aircraft with flaps)
    % CL_min = 0.5;             % Minimum CL (reasonable lower bound)
    % 
    % % Iterate until the velocity is trimmed to the target
    % iteration = 0;
    % while abs(V_current - target_velocity) > tolerance
    %     % Calculate dynamic lift coefficient (CL) based on AoA and flap settings
    %     % Here AoA is replaced by glide path angle (gamma), which will control the desired CL
    % 
    %     % Calculate the current lift (L) and drag (D) forces
    %     alpha = atan(W / (0.5 * rho * V_current^2 * S * C_L0));  % Initial alpha calculation
    % 
    %     % Calculate CL using current angle of attack
    %     C_L = C_L0 + C_L_alpha * alpha;
    % 
    %     % Ensure CL stays within reasonable limits (0.5 < CL < 2.5)
    %     if C_L > CL_max
    %         C_L = CL_max;
    %     elseif C_L < CL_min
    %         C_L = CL_min;
    %     end
    % 
    %     % Induced drag model (quadratic relation to lift coefficient)
    %     C_D = C_D0 + k * C_L^2;
    % 
    %     % Calculate steady glide conditions
    %     L = 0.5 * rho * V_current^2 * C_L * S;  % Lift force (N)
    %     D = 0.5 * rho * V_current^2 * C_D * S;  % Drag force (N)
    % 
    %     % Calculate glide path angle (gamma) for current L and D
    %     gamma = atan(D / L);  % Glide path angle (radians)
    % 
    %     % Angle of attack is approximately equal to the glide path angle
    %     alpha = gamma;
    % 
    %     % Calculate distance traveled during this step
    %     distance_traveled = distance_traveled + V_current * 1;  % Assume each iteration is 1 second
    % 
    %     % Apply drag force to reduce velocity
    %     V_new = V_current - D / m * 1;  % Apply drag force to reduce speed (assume 1 second steps)
    % 
    %     % Update current velocity and ensure it stays above the target
    %     V_current = max(V_new, target_velocity);
    % 
    %     % Prevent angle of attack from reaching stall (ensure alpha does not exceed stall angle)
    %     if alpha > stall_angle
    %         alpha = stall_angle;
    %     end
    % 
    %     % Increment iteration count
    %     iteration = iteration + 1;
    % 
    %     % Display progress
    %     if mod(iteration, 10) == 0
    %         fprintf('Iteration %d: Velocity = %.2f m/s, Distance Traveled = %.2f m, CL = %.2f, Alpha = %.2f deg\n', ...
    %             iteration, V_current, distance_traveled, C_L, rad2deg(alpha));
    %     end
    % end
    % 
    % % Final Output
    % fprintf('Final Trimming Results:\n');
    % fprintf('Final Glide Velocity: %.2f m/s\n', V_current);
    % fprintf('Total Distance Traveled: %.2f m\n', distance_traveled);
    % fprintf('Final Lift Coefficient (CL): %.2f\n', C_L);
    % fprintf('Final Angle of Attack (Alpha): %.2f degrees\n', rad2deg(alpha));
    % 
    % % Plotting the velocity reduction over time
    % figure;
    % plot(1:iteration, linspace(100, V_current, iteration), 'LineWidth', 2);
    % xlabel('Iteration');
    % ylabel('Velocity (m/s)');
    % title('Reduction in Glide Velocity Over Time');
    % grid on;


    % Aircraft parameters
    m = 5000;                % Mass of the aircraft (kg)
    S = 30;                  % Wing reference area (m²)
    C_L0 = 0.3;              % Lift coefficient at zero angle of attack (dimensionless)
    C_L_alpha = 0.1;         % Lift curve slope (per radian) (dimensionless)
    C_D0 = 0.02;             % Parasitic drag coefficient (dimensionless)
    k = 0.04;                % Induced drag factor (dimensionless)
    g = 9.81;                % Gravitational acceleration (m/s²)
    
    % Initial conditions
    altitude = 3000;         % Initial altitude (m)
    V_current = 100;         % Initial horizontal speed (m/s)
    target_velocity = 70;    % Target glide velocity (m/s) - optimize this value for range
    time_step = 1;           % Time step for simulation (seconds)
    distance_traveled = 0;   % Distance traveled (m)
    total_velocity = V_current;  % Total velocity at start
    glide_angle = 0;         % Initial glide angle (degrees)
    
    % Initialize arrays to store results
    altitudes = altitude;    % Altitude over time
    distances = distance_traveled; % Horizontal distance over time
    velocities = total_velocity; % Speed over time
    glide_angles = glide_angle;  % Glide angle over time

    % Iterate until the velocity gets close to the target value
    while V_current > target_velocity
        % Calculate air density at current altitude using the ISA model
        [~, rho, ~, ~] = atmosisa(altitude);  % Air density (kg/m³)
        
        % Calculate the dynamic pressure
        q = 0.5 * rho * V_current^2;  % Dynamic pressure
        
        % Calculate lift coefficient based on angle of attack (alpha)
        alpha = atan(g / V_current);  % Assuming gliding with minimal angle of attack
        C_L = C_L0 + C_L_alpha * alpha;
        
        % Induced drag coefficient (function of CL)
        C_D = C_D0 + k * C_L^2;
        
        % Lift and Drag forces
        L = q * C_L * S;  % Lift force (N)
        D = q * C_D * S;  % Drag force (N)
        
        % Update glide angle (gamma) based on lift and drag
        glide_angle = atan(D / L);  % Glide path angle (radians)
        
        % Update altitude and distance traveled
        altitude = altitude - V_current * sin(glide_angle) * time_step;  % Altitude decreases as aircraft glides down
        distance_traveled = distance_traveled + V_current * cos(glide_angle) * time_step;  % Horizontal distance
        
        % Update velocity (bleeding speed due to drag)
        V_current = V_current - D / m * time_step;  % Update velocity using drag force
        
        % Store values for plotting
        altitudes = [altitudes, altitude];
        distances = [distances, distance_traveled];
        velocities = [velocities, V_current];
        glide_angles = [glide_angles, rad2deg(glide_angle)];
    end
    
    % Final Outputs
    fprintf('Final Glide Simulation Results:\n');
    fprintf('Final Glide Velocity: %.2f m/s\n', V_current);
    fprintf('Total Distance Traveled: %.2f m\n', distance_traveled);
    fprintf('Final Altitude: %.2f m\n', altitude);
    fprintf('Final Glide Angle: %.2f degrees\n', rad2deg(glide_angle));
    
    % Plotting Results
    figure;
    subplot(3,1,1);
    plot(distances, altitudes);
    xlabel('Distance (m)');
    ylabel('Altitude (m)');
    title('Altitude vs. Horizontal Distance');

    subplot(3,1,2);
    plot(distances, velocities);
    xlabel('Distance (m)');
    ylabel('Velocity (m/s)');
    title('Velocity vs. Horizontal Distance');
    
    subplot(3,1,3);
    plot(distances, glide_angles);
    xlabel('Distance (m)');
    ylabel('Glide Angle (degrees)');
    title('Glide Angle vs. Horizontal Distance');

