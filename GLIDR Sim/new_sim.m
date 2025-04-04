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


config = 2; % Change to run sims on each config
% 1 ==== Blended Wing Body
% 2 ==== Rogallo Wing
% 3 ==== Bullet Bill


DESIGN = configuration(config);

end_sim = size(DESIGN.S,1);
for ii = 1:1:end_sim



% Freefall until desired q
options = odeset('Events', @y1_free);
place_free = @(t,y) funfree(t,y,DESIGN,ii);
[t,y,te,ye,ei] = ode45(place_free, linspace(0,100,1000), [36000;0], options); % need to impliment aoa for linear region

free(:,:,ii) = {t;y(:,1);abs(y(:,2))};
free_out(:,:,ii) = [te,ye,ei];

global pull_init
pull_init = te;

pull_inip(1,ii) = te;

% Pullout @ desired q
options = odeset('Events', @y1_pull);
place_pull = @(t,y) funpull(t,y,DESIGN,ii);
[t,y,te,ye,ei] = ode45(place_pull, [0 5000], [0;ye(1,1);ye(1,2);-pi/2], options); % need to impliment theta and Q

pull(:,:,ii) = {t;y(:,1);y(:,2);abs(y(:,3));y(:,4)};
pull_out(:,:,ii) = [te,ye,ei];

glide_initp(1,ii) = te;

% going from fast to ideal cruise

% options = odeset('Events', @y1_glide);
% chat = @(t,y) funglide(t,y,DESIGN,ii);
% [t,y,te,ye,ei] = ode15s(chat, [0 100], [ye(1,1);ye(1,2);abs(ye(1,3)); 0], options); % need to impliment theta and Q

ic = [ye(1,1);ye(1,2);abs(ye(1,3)); 0];
[t,y,L_D,ic,te] = sadglide(ic,DESIGN,ii);

%velo = cell(size(t,1),2);
velo(:,1) = sqrt(y(:,3).^2+y(:,4).^2);
%t = zeros(size(velo,1),1);
glide(:,:,ii) = {t;y(:,1);y(:,2);velo(:,1)};

clear velo
% Stead state cruise

[t,y,L_D] = funsteady(ic,DESIGN,ii);

velo(:,1) = sqrt(y(:,3).^2+y(:,4).^2);
steady(:,:,ii) = {t;y(:,1);y(:,2);velo(:,1)};
steady_inip(1,ii) = te; 

clear velo



% CL DESIGN POINT EQUATION ++++++++++++++++++++++++++++++++++++++++++++++++

% cL_glide = sqrt(3.*cD0.*(pi.*DESIGN.e.*DESIGN.AR));
% V_trim = sqrt(2.*DESIGN.m.*DESIGN.g./(rho.*DESIGN.S(ii)).*sqrt((1/(pi.*DESIGN.e(ii).*DESIGN.AR(ii)))/(3.*CD0)));
% gamma = atan(CD/CL_glide);
% Glide for min sink^^

% Landing (Bleed V for safe landing)

end


%% Cool Plots

 i = 1;
for free_ind = 1:size(free,3)
    figure(i);
    free_mat_ind = cell2mat(free(:,:,free_ind));
    t = free_mat_ind(1:size(free_mat_ind,1)/3);
    alt = free_mat_ind(size(free_mat_ind,1)/3+1:size(free_mat_ind,1)/3*2);
    speed = free_mat_ind(size(free_mat_ind,1)/3*2+1:end);
    hold on;
    plot(t,alt);
    title("alt vs time")
    legend();
    i = i + 1;
    figure(i);
    hold on;
    plot(speed,alt);
    title("alt vs speed")
    legend();
    i = i - 1;
end
for pull_ind = 1:size(pull,3)
    figure(i);
    pull_mat_ind = cell2mat(pull(:,:,pull_ind));
    t = pull_mat_ind(1:size(pull_mat_ind,1)/5);
    x = pull_mat_ind(size(pull_mat_ind,1)/5+1:size(pull_mat_ind,1)/5*2);
    alt = pull_mat_ind(size(pull_mat_ind,1)/5*2+1:size(pull_mat_ind,1)/5*3);
    speed = pull_mat_ind(size(pull_mat_ind,1)/5*3+1:size(pull_mat_ind,1)/5*4);
    hold on;
    plot(t+pull_inip(pull_ind),alt);
    title("alt vs time")
    legend();
    i = i + 1;
    figure(i);
    hold on;
    plot(speed,alt);
    title("alt vs speed")
    legend();
    i = i + 1;
    figure(i);
    hold on;
    plot(x,alt);
    title("alt vs x-dist")
    legend();
    i = i - 2;
end
for glide_ind = 1:size(glide,3)
    figure(i);
    glide_mat_ind = cell2mat(glide(:,:,glide_ind));
    t = glide_mat_ind(1:size(glide_mat_ind,1)/4);
    x = glide_mat_ind(size(glide_mat_ind,1)/4+1:size(glide_mat_ind,1)/4*2);
    alt = glide_mat_ind(size(glide_mat_ind,1)/4*2+1:size(glide_mat_ind,1)/4*3);
    speed = glide_mat_ind(size(glide_mat_ind,1)/4*3+1:end);
    hold on;
    plot(t+glide_initp(glide_ind)+pull_inip(glide_ind),alt);
    title("alt vs time")
    legend();
    i = i + 1;
    figure(i);
    hold on;
    plot(speed,alt);
    title("alt vs speed")
    legend();
    i = i + 1;
    figure(i);
    hold on;
    plot(x,alt);
    title("alt vs x-dist")
    legend();
    i = i - 2;
end
for steady_ind = 1:size(steady,3)
    figure(i);
    steady_mat_ind = cell2mat(steady(:,:,steady_ind));
    t = steady_mat_ind(1:size(steady_mat_ind,1)/4);
    x = steady_mat_ind(size(steady_mat_ind,1)/4+1:size(steady_mat_ind,1)/4*2);
    alt = steady_mat_ind(size(steady_mat_ind,1)/4*2+1:size(steady_mat_ind,1)/4*3);
    speed = steady_mat_ind(size(steady_mat_ind,1)/4*3+1:end);
    hold on;
    plot(t+glide_initp(steady_ind)+pull_inip(steady_ind)+steady_inip(steady_ind),alt);
    title("alt vs time")
    legend();
    i = i + 1;
    figure(i);
    hold on;
    plot(speed,alt);
    title("alt vs speed")
    legend();
    i = i + 1;
    figure(i);
    hold on;
    plot(x,alt);
    title("alt vs x-dist")
    legend();
    i = i - 2;
end




% Plotting tools
% plot(free(:,1,:),free(:,1,:))



%% Functions for simulations


% Freefall section ode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function dydt = funfree(t,y,DESIGN,ii)


% Get atmospheric properties at current height
    [T, a, P, rho, nu, mu] = atmosisa(y(1), 'extended', true);

    % Drag Buildup for changing Re and flight condition
    Re = rho.*DESIGN.c.*y(2)./mu;


cL = 0;
[cD0,cDi] = drag(Re,cL,DESIGN,ii);
cD = cD0 + cDi; % Add flap thing drag

drag_force = 0.5 * cD * DESIGN.S(ii) * rho * y(2).^2;

    % Calculate the drag force and acceleration due to gravity
     % Drag force (assumes velocity is y(2))
    dydt = [y(2); -DESIGN.g + drag_force ./ DESIGN.m]; % [dy/dt, dv/dt]
end

function [position,isterminal,direction] = y1_free(t,y)

position = y(2) + 200; % The value that we want to be zero (speed)
    isterminal = 1;  % Halt integration 
    direction = 0;   % The zero can be approached from either direction
end



% Pullout section ode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function dydt = funpull(t,y,DESIGN,ii)



% Get atmospheric properties at current height
    [T, a, P, rho, nu, mu] = atmosisa(y(2), 'extended', true);

    % Drag Buildup for changing Re and flight condition
    Re = rho.*DESIGN.c.*abs(y(3))./mu;

    % TO BE FURTHER EVALUATED++++++++++++++++++++++++++++++++++++++++
    % Estimation based on flaps to to create negative lift 
    global pull_init

    if t == pull_init
        cL = 0;
    else
        cL = -0.5;
    end
   
    
[cD0,cDi] = drag(Re,cL,DESIGN,ii);
cD = cD0 + cDi;
    
drag_force = 0.5 .* rho .* y(3)^2 .* DESIGN.S(ii) .* cD;
lift_force = 0.5 .* rho .* y(3)^2 .* DESIGN.S(ii) .* cL;
   
    
% Q assuming aircraft is point mass. Can be updated with inertia. Replace m
% with Izz to do so. ++++++++++++++++++++++++++++++++++++++++++++++++++++

% Assumes that aero center is on same waterline as cg.++++++++++++++++++++

Q = ( .5*rho*(y(3)^2)*cL*(DESIGN.S(ii)/DESIGN.m) - DESIGN.g*cos(y(4)) )/y(3);
    
 % Calculate the drag force and acceleration due to gravity
     % Drag force (assumes velocity is y(2))
    dydt = [abs(y(3)).*cos(y(4)); abs(y(3)).*sin(y(4)); DESIGN.g.*sin(y(4))-drag_force/DESIGN.m; Q]; % [dx/dt, dy/dt, dv/dt, dtheta/dt]
end

function [position,isterminal,direction] = y1_pull(t,y)
    position = y(4); % The value that we want to be zero (theta)
    isterminal = 1;  % Halt integration 
    direction = 0;   % The zero can be approached from either direction
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% function [t,y,L_D] = sadglide(ic,DESIGN,ii)
% 
% 
% dh = 1000;
% ind = round(ic(2)/dh);
% y = zeros(ind,4);
% 
% for jj  = 1:ind
%     [T, a, P, rho, nu, mu] = atmosisa(y(jj,2), 'extended', true);
% 
%     if jj == 1
%         V = sqrt(ic(3)^2+ic(4)^2);
%     else
%         V = sqrt(y(jj-1,3)^2+y(jj-1,4)^2);
%     end
%     W = DESIGN.m.*DESIGN.g;
%     gamma = -4.*pi/180;
% 
%     Re = rho.*DESIGN.c.*V./mu;
% 
% y(jj,5) = (rho.*V^2.*DESIGN.S(ii));
% cL = 2.*W.*cos(gamma)./y(jj,5);
% [cD0,cDi] = drag(Re,cL,DESIGN,ii);
% cD = cD0 + cDi;
% 
% L_D(1,jj) = cL/cD;
% dx(1,jj) = dh.*L_D(1,jj);
% 
% if jj == 1
%     y(jj,1) = dx(1,jj);
%     y(jj,2) = ic(2);
% else
% y(jj,1) = y(jj-1,1) + dx(1,jj);
% y(jj,2) = ic(2) - jj*dh;
% end
% 
% t(jj,1) = V.*sqrt(dx(1,jj)^2+dh^2);
% 
% if jj == 1
%     y(jj,3) = ic(3);
%     y(jj,4) = ic(4);
% else
% y(jj,3) = y(jj-1,3) / dx(1,jj);
% y(jj,4) = y(jj-1,4) / dh;
% end
% 
% % V(jj,1) = 
% 
% 
% end
% end

% function [t,y,L_D] = sadglide(ic,DESIGN,ii)
%     %dh = 1000;
%     %ind = round(ic(2)/dh);
%     y = zeros(1000,4);  % Initialize result matrix
% 
%     % Define time step (dt) for velocity update
%     dt = 0.1;  % Time step (adjust if needed)
% 
%     alt = ic(1);
% 
%     for jj = 1:1000
%         [T, a, P, rho, nu, mu] = atmosisa(alt, 'extended', true);
% 
%         % Initial velocity
%         if jj == 1
%             V = sqrt(ic(3)^2 + ic(4)^2);
%         else
%             V = sqrt(y(jj-1,3)^2 + y(jj-1,4)^2);
%         end
%         dh = y(ii,4).*dt;
%         W = DESIGN.m * DESIGN.g;
%         gamma = -4 * pi / 180;  % Glide angle
% 
%         Re = rho * DESIGN.c * V / mu;  % Reynolds number
% 
%         denom = (rho * V^2 * DESIGN.S(ii));  % Dynamic pressure * S
%         cL = 2 * W * cos(gamma) / denom;  % Lift coefficient
%         [cD0, cDi] = drag(Re, cL, DESIGN, ii);  % Drag coefficients
%         cD = cD0 + cDi;  % Total drag coefficient
% 
%         L_D(1,jj) = cL / cD;
%         dx(1,jj) = dh * L_D(1,jj);
% 
%         if jj == 1
%             y(jj,1) = dx(1,jj);  % Horizontal position
%             y(jj,2) = ic(2);  % Vertical position
%         else
%             y(jj,1) = y(jj-1,1) + dx(1,jj);  % Update horizontal position
%             y(jj,2) = y(jj-1,2) - dh;  % Update vertical position
%         end
% 
%         alt = y(jj,2);
%         % Update velocity based on aerodynamic forces
%         D = 0.5 * rho * V^2 * DESIGN.S(ii) * cD;  % Drag force
%         L = 0.5 * rho * V^2 * DESIGN.S(ii) * cL;  % Lift force
% 
%         % Accelerations
%         ax = -D * cos(gamma) / DESIGN.m;  % Horizontal acceleration
%         ay = -W + L * sin(gamma) / DESIGN.m;  % Vertical acceleration
% 
%         % Update velocities
%         if jj == 1
%             y(jj,3) = ic(3);  % Initial horizontal velocity
%             y(jj,4) = ic(4);  % Initial vertical velocity
%         else
%             y(jj,3) = y(jj-1,3) + ax * dt;  % Horizontal velocity update
%             y(jj,4) = y(jj-1,4) + ay * dt;  % Vertical velocity update
%         end
% 
%         % Time step (optional for display)
%         t(jj,1) = V * sqrt(dx(1,jj)^2 + dh^2);  % Time step calculation (you can adjust this if needed)
%     end
% end

function [t,y,L_D,ic,te] = sadglide(ic,DESIGN,ii)
% Aircraft parameters

    C_L0 = 0.3;              % Lift coefficient at zero angle of attack (dimensionless)
             % Parasitic drag coefficient (dimensionless)
    k = 1/(DESIGN.e(ii)*pi*DESIGN.AR(ii));                % Induced drag factor (dimensionless)

    g = 9.81;                % Gravitational acceleration (m/s²)
    
    % Flight conditions (inputs)
    alt = ic(2);         % Altitude (m)
    [T, a, P, rho, nu, mu] = atmosisa(alt, 'extended', true);

    % Drag Buildup for changing Re and flight condition
    Re = rho.*DESIGN.c.*abs(ic(3))./mu;
    cD0 = drag(Re,0,DESIGN,ii);

    % Compute weight
    W = DESIGN.m * DESIGN.g;               % Weight of the aircraft (N)
    
    % Stall angle (radians) (stall angle assumed for typical subsonic aircraft)
    stall_angle = deg2rad(15);  
    
    % Target glide velocity (optimal velocity determined earlier)
    % Assuming the optimized value was calculated previously
    target_velocity = sqrt((2 * W) / (rho * DESIGN.S(ii) * C_L0));     % Example target velocity in m/s (can be calculated from previous results)
    
    % Initial condition: Aircraft starts with a higher glide velocity
    V_current = abs(ic(3));           % Higher initial glide velocity (m/s)
    tolerance = 0.01;         % Convergence tolerance for velocity
    xdistance_traveled = 0;    % Distance traveled (m)
    ydistance_traveled = 0;
    
    % Define lift coefficient limits (realistic values)
    CL_max = 2.5;             % Maximum CL (typically around this for subsonic aircraft with flaps)
    CL_min = 0.1;             % Minimum CL (reasonable lower bound)

    % Iterate until the velocity is trimmed to the target
    iter = 0;
    while abs(V_current - target_velocity) > tolerance
        % Calculate dynamic lift coefficient (CL) based on AoA and flap settings
        % Here AoA is replaced by glide path angle (gamma), which will control the desired CL
        
        % Calculate the current lift (L) and drag (D) forces
        alpha = atan(W / (0.5 * rho * V_current^2 * DESIGN.S(ii) * C_L0));  % Initial alpha calculation
        
        % Calculate CL using current angle of attack
        % C_L = C_L0 + DESIGN.cLalpha(ii) * alpha
        C_L = W / (0.5 * rho * V_current^2 * DESIGN.S(ii));
        
        % Ensure CL stays within reasonable limits (0.5 < CL < 2.5)
        if C_L > CL_max
            C_L = CL_max;
        elseif C_L < CL_min
            C_L = CL_min;
        end
        
        % Induced drag model (quadratic relation to lift coefficient)
        C_D = cD0 + k * C_L^2;
        
        % Calculate steady glide conditions
        L = 0.5 * rho * V_current^2 * C_L * DESIGN.S(ii);  % Lift force (N)
        D = 0.5 * rho * V_current^2 * C_D * DESIGN.S(ii) * DESIGN.S_ratio;  % Drag force (N)
        
% INSTEAD BLEEDING SPEED AT A CONSTANT ALT================================

        % Calculate glide path angle (gamma) for current L and D
        gamma = D / L;  % Glide path angle (radians)
        
        % Angle of attack is approximately equal to the glide path angle
        % alpha = gamma;
        
        % Velocity slip into x and y components for distance
        V_Cx = V_current.*cos(gamma);
        V_Cy = V_current.*sin(gamma);
        V_Cy = 0;
        % Calculate distance traveled and altitude lost during this step
        xdistance_traveled = + V_Cx * 1;  % Assume each iteration is 1 second
        ydistance_traveled = + V_Cy * 1;  % Assume each iteration is 1 second

        % Apply drag force to reduce velocity
        V_new = V_current - D / DESIGN.m * 1;  % Apply drag force to reduce speed (assume 1 second steps)
        
        % Update current velocity and ensure it stays above the target
        V_current = max(V_new, target_velocity);
        
        % Prevent angle of attack from reaching stall (ensure alpha does not exceed stall angle)
        if alpha > stall_angle
            alpha = stall_angle;
        end
        
        % Increment iteration count and saving needed parameters
        iter = iter + 1;
        if iter == 1
        L_D(iter,1) = L/D; 
        y(iter,1) = ic(1);
        y(iter,2) = alt;
        y(iter,3) = ic(3);
        y(iter,4) = V_Cy;
        t(iter,1) = iter;
        else
        L_D(iter,1) = L/D; 
        y(iter,1) = y(iter-1,1) + xdistance_traveled;
        y(iter,2) = y(iter-1,2) - ydistance_traveled;
        y(iter,3) = V_Cx;
        y(iter,4) = V_Cy;
        t(iter,1) = iter;
        
    end
    end

    ic = [y(iter,1), y(iter,2), y(iter,3), y(iter,4)];
    te = t(iter);
end


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% function dydt = funglide(t,y,DESIGN,ii)
% 
% % Get atmospheric properties at current height
%     [T, a, P, rho, nu, mu] = atmosisa(y(1), 'extended', true);
% 
%     V = sqrt(y(3)^2+y(4)^2)
%     % Drag Buildup for changing Re and flight condition
%     Re = rho.*DESIGN.c.*V./mu;
% 
% 
%     W = DESIGN.m.*DESIGN.g;
% 
%     gamma = -4.*pi/180;
% 
% cL = 2.*W.*cos(gamma)./(rho.*DESIGN.S(ii).*V^2);
% [cD0W,cDiW] = drag(Re,cL,DESIGN,ii);
% 
% cLH = DESIGN.S(ii).*cL.*DESIGN.M_arm/(DESIGN.SH.*DESIGN.M_armH);
% [cD0H,cDiH] = drag(Re,cLH,DESIGN,ii);
% 
% 
% % V_trim_g = sqrt(2.*DESIGN.m.*DESIGN.g./(rho.*DESIGN.S(ii)).*sqrt((1/(pi.*DESIGN.e(ii).*DESIGN.AR(ii)))/(3.*cD0)));
% % cLW = sqrt(3.*cD0.*(pi.*DESIGN.e.*DESIGN.AR));
% 
% % x_velo = velo.*cos(gamma);
% % y_velo = velo.*sin(gamma);
% 
% 
% cD = cD0W + cDiW + cD0H + cDiH
% drag_force = 0.5 .* rho .* V^2 .* DESIGN.S(ii) .* cD;
% lift_force = 0.5 .* rho .* V^2 .* DESIGN.S(ii) .* cL;
% %V_trim_g = sqrt(2.*DESIGN.m.*DESIGN.g./(rho.*DESIGN.S(ii)).*sqrt((1/(pi.*DESIGN.e(ii).*DESIGN.AR(ii)))/(3.*cD0)));
% 
% 
% 
% x_accel = -drag_force.*cos(gamma)/DESIGN.m - lift_force.*sin(gamma)/DESIGN.m;
% y_accel =  -drag_force.*sin(gamma)/DESIGN.m + lift_force.*cos(gamma)/DESIGN.m - DESIGN.g;
% 
% dydt = [y(3); y(4); -drag_force.*cos(gamma)/DESIGN.m - lift_force.*sin(gamma)/DESIGN.m; -drag_force.*sin(gamma)/DESIGN.m + lift_force.*cos(gamma)/DESIGN.m - DESIGN.g];
% % Glide for min sink^^
% 
% end


function [position,isterminal,direction] = y1_glide(t,y)
    position = y(2); % The value that we want to be zero (altitude)
    isterminal = 1;  % Halt integration 
    direction = 0;   % The zero can be approached from either direction
end

function [t,y,L_D] = funsteady(ic,DESIGN,ii)
% Glide has constant q to keep aero forces
alt = ic(2);

iter = 0;
y(1,1) = ic(1);
y(1,2) = ic(2);
y(1,3) = ic(3);
y(1,4) = ic(4);

while alt > 1700
    iter = iter + 1; % One iteration is 5 second

     [T, a, P, rho, nu, mu] = atmosisa(alt, 'extended', true);

     if iter == 1
     % finding set q value
     q_con = 0.5*rho*y(iter,3)^2;
     end
    
     V = sqrt(2*q_con/(rho));
     W = DESIGN.m*DESIGN.g;
     
    cL = W / (0.5 * rho * V^2 * DESIGN.S(ii));
   

    % Drag Buildup for changing Re and flight condition
    Re = rho.*DESIGN.c.*V./mu;
    [cD0,cDi] = drag(Re,cL,DESIGN,ii);
    cD = (cD0 + cDi);

    L_Dtemp = cL/cD;
    gamma = atan(cD/cL);

    if iter == 1
        L_D(iter,1) = L_Dtemp;
    else
        y(iter,3) = V*cos(gamma);
        y(iter,4) = V*sin(gamma);
        y(iter,1) = y(iter-1,1) + y(iter,3) * 5;
        y(iter,2) = y(iter-1,2) - y(iter,4) * 5;
        t(iter,1) = iter * 5;
        L_D(iter,1) = L_Dtemp;
    end
    alt = y(iter,2);
end


end