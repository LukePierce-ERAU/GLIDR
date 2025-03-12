%% Full simulation of fall from 36.5km for DCR configurations
clear
close all
%% GEOMETRY ASSUMPTION SECTIONS AND GLOBAL VARIABLES ++++++++++++++++++++++


config = 3; % Change to run sims on each config
% 1 ==== Blended Wing Body
% 2 ==== Rogallo Wing
% 3 ==== Bullet Bill


DESIGN = configuration(config);
steady = struct();

end_sim = size(DESIGN.S,2);
for ii = 1:1:end_sim

% Freefall until desired q
options = odeset('Events', @y1_free);
place_free = @(t,y) funfree(t,y,DESIGN,ii);
[t,y,te,ye,ei] = ode45(place_free, [0 10000], [36000;0], options); % need to impliment aoa for linear region

for k = 1:numel(t)
    [~,moment(k,:)] = place_free(t(k),y(k,:));
end

freevar(:,:,ii) = {t;y(:,1);abs(y(:,2))};
free_out(:,:,ii) = [te,ye,ei];

global pull_init
pull_init = te;

glide_inip(1,ii) = te;
ic = [0;ye(1,1);abs(ye(1,2)); 0];
clear velo

% Stead state cruise
[t,y,L_D,cL] = funsteady(ic,DESIGN,ii);

velo(:,1) = sqrt(y(:,3).^2+y(:,4).^2);
steadyvar(:,:,ii) = {t;y(:,1);y(:,2);velo(:,1);L_D(:,1);cL(:,1)};
steady_inip(1,ii) = te; 

m_name = ['moment', num2str(ii)];
steady.(m_name) = moment;

clear velo
end


%% Cool Plots

free = struct();
master = struct();

for free_ind = 1:size(freevar,3)
    free_mat_ind = cell2mat(freevar(:,:,free_ind));
    t = free_mat_ind(1:size(free_mat_ind,1)/3);
    x = zeros(size(free_mat_ind,1)/3,1);
    alt = free_mat_ind(size(free_mat_ind,1)/3+1:size(free_mat_ind,1)/3*2);
    speed = free_mat_ind(size(free_mat_ind,1)/3*2+1:end);
    t_name = ['t', num2str(free_ind)];
    x_name = ['x', num2str(free_ind)];
    alt_name = ['alt', num2str(free_ind)];
    speed_name = ['speed', num2str(free_ind)];
    free.(t_name) = t;
    free.(x_name) = x;
    free.(alt_name) = alt;
    free.(speed_name) = speed;
end



for steady_ind = 1:size(steadyvar,3)
    steady_mat_ind = cell2mat(steadyvar(:,:,steady_ind));
    t = steady_inip(1,steady_ind) + steady_mat_ind(1:size(steady_mat_ind,1)/6);
    x = steady_mat_ind(size(steady_mat_ind,1)/6+1:size(steady_mat_ind,1)/6*2);
    alt = steady_mat_ind(size(steady_mat_ind,1)/6*2+1:size(steady_mat_ind,1)/6*3);
    speed = steady_mat_ind(size(steady_mat_ind,1)/6*3+1:size(steady_mat_ind,1)/6*4);
    L_D = steady_mat_ind(size(steady_mat_ind,1)/6*4+1:size(steady_mat_ind,1)/6*5);
    cL = steady_mat_ind(size(steady_mat_ind,1)/6*5+1:end);
    % Add pulls for L/D ratio and cL
    t_name = ['t', num2str(steady_ind)];
    x_name = ['x', num2str(steady_ind)];
    alt_name = ['alt', num2str(steady_ind)];
    speed_name = ['speed', num2str(steady_ind)];
    L_D_name = ['LoveD', num2str(steady_ind)];
    cL_name = ['cL', num2str(steady_ind)];
    steady.(t_name) = t;
    steady.(x_name) = x;
    steady.(alt_name) = alt;
    steady.(speed_name) = speed;
    steady.(L_D_name) = L_D;
    steady.(cL_name) = cL;
end

% Assign S and AR to Each Flight Configuration
configs = {'free', 'pull', 'glide', 'steady'};

for k = 1:length(configs)
    for jj = 1:length(DESIGN.S)  % Ensure indexing matches DESIGN.S
        eval([configs{k}, '.([''S'', num2str(jj)]) = DESIGN.S(jj);']);  % Assign S values
        eval([configs{k}, '.([''AR'', num2str(jj)]) = DESIGN.AR(jj);']);  % Assign AR values
    end
end

% Creating one line for each config

for jj = 1:1:ii
master.(['t', num2str(jj)]) = vertcat(free.(['t',num2str(jj)]),steady.(['t',num2str(jj)]));
master.(['x', num2str(jj)]) = vertcat(free.(['x',num2str(jj)]),steady.(['x',num2str(jj)]));
master.(['alt', num2str(jj)]) = vertcat(free.(['alt',num2str(jj)]),steady.(['alt',num2str(jj)]));
master.(['speed', num2str(jj)]) = vertcat(free.(['speed',num2str(jj)]),steady.(['speed',num2str(jj)]));
master.(['LoveD', num2str(jj)]) = vertcat(steady.(['LoveD',num2str(jj)]));
master.(['cL', num2str(jj)])    = vertcat(steady.(['cL',num2str(jj)]));

% Add S and AR using vertcat
master.(['S', num2str(jj)]) = vertcat(free.(['S', num2str(jj)]), pull.(['S', num2str(jj)]), glide.(['S', num2str(jj)]), steady.(['S', num2str(jj)]));
master.(['AR', num2str(jj)]) = vertcat(free.(['AR', num2str(jj)]), pull.(['AR', num2str(jj)]), glide.(['AR', num2str(jj)]), steady.(['AR', num2str(jj)]));

end


%% Range vs. Wing Area (S)
figure
hold on
plot([master.S1, master.S2, master.S3, master.S4, master.S5, ...
      master.S6, master.S7, master.S8, master.S9, master.S10], ...
     [master.x1(end)/1000, master.x2(end)/1000, master.x3(end)/1000, master.x4(end)/1000, master.x5(end)/1000, ...
      master.x6(end)/1000, master.x7(end)/1000, master.x8(end)/1000, master.x9(end)/1000, master.x10(end)/1000], '-o', 'LineWidth', 2)
grid on
title('Range vs. Wing Area', 'FontSize', 18)
xlabel('Wing Area (S) [m²]', 'FontSize', 16)
ylabel('Range [km]', 'FontSize', 16)

%% Range vs. Aspect Ratio (AR)
figure
hold on
plot([master.AR1, master.AR2, master.AR3, master.AR4, master.AR5, ...
      master.AR6, master.AR7, master.AR8, master.AR9, master.AR10], ...
     [master.x1(end)/1000, master.x2(end)/1000, master.x3(end)/1000, master.x4(end)/1000, master.x5(end)/1000, ...
      master.x6(end)/1000, master.x7(end)/1000, master.x8(end)/1000, master.x9(end)/1000, master.x10(end)/1000], '-o', 'LineWidth', 2)
grid on
title('Range vs. Aspect Ratio', 'FontSize', 18)
xlabel('Aspect Ratio (AR)', 'FontSize', 16)
ylabel('Range [km]', 'FontSize', 16)

%% Altitude vs Lift to Drag ratio
figure
hold on
plot(master.cL1(2:end),vertcat(steady.alt1(2:end))/1000,'b','LineWidth',2)
plot(master.cL6(2:end),vertcat(steady.alt6(2:end))/1000,'k','LineWidth',2)
plot(master.cL10(2:end),vertcat(steady.alt10(2:end))/1000,'m','LineWidth',2)
grid on
title('Altitude over Target c_{L}', 'FontSize',18)
xlabel('c_{L}', 'FontSize',16)
ylabel('Altitude [km]', 'FontSize',16)

%% Altitude over Time
figure
hold on
plot(master.t1/60,master.alt1/1000,'m','LineWidth',2)
plot(master.t6/60,master.alt6/1000,'k','LineWidth',2)
plot(master.t10/60,master.alt10/1000,'b','LineWidth',2)
grid on
title('Altitude over Time','FontSize',18)
xlabel('Time [mins]','FontSize',16)
ylabel('Altitude [km]','FontSize',16)

%% 2D Flight Path
figure
hold on
plot(master.x1/1000,master.alt1/1000,'m','LineWidth',2)
plot(master.x6/1000,master.alt6/1000,'k','LineWidth',2)
plot(master.x10/1000,master.alt10/1000,'b','LineWidth',2)
xline(45,'r')
grid on
title('2D Flight path','FontSize',18)
xlabel('Distance Travelled [km]','FontSize',16)
ylabel('Altitude [km]','FontSize',16)

%% Altitude VS Speed
figure
hold on
plot(master.speed1,master.alt1/1000,'m','LineWidth',2)
plot(master.speed6,master.alt6/1000,'k','LineWidth',2)
plot(master.speed10,master.alt10/1000,'b','LineWidth',2)
grid on
title('Altitude over Speed','FontSize',18)
xlabel('Speed [m/s]','FontSize',16)
ylabel('Altitude [km]','FontSize',16)

%% Figure: Airbreak moment vs. Altitude
figure;
plot(steady.moment1, free.alt1, 'b', 'LineWidth', 2);
set(gca, 'YDir', 'reverse'); % Altitude decreases downward
xlabel('Moment');
ylabel('Altitude (m)');
title('Moment vs. Altitude');
grid on;

%% Functions for simulations

% Freefall section ode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [dydt,moment] = funfree(t,y,DESIGN,ii)


% Get atmospheric properties at current height
    [T, a, P, rho, nu, mu] = atmosisa(y(1), 'extended', true);

    q = 0.5 * rho * y(2)^2; % Dynamic pressure
    
    % Drag Buildup for changing Re and flight condition
    Re = rho.*DESIGN.c.*abs(y(2))./mu;
    cL = 0;
    [cD0,cDi] = drag(Re,cL,DESIGN,ii);
    cD = cD0 + cDi;

    C_d_airbrake = 2.2; % Given Cd
    D_airbrake = q *  DESIGN.S_airbrake * C_d_airbrake; % Airbrake drag force
    moment(numel(t),1) = D_airbrake * DESIGN.airbrake_ac; % This is for 1 of 3 airbrakes

    drag_force = 0.5 * cD * DESIGN.S(ii) * rho * y(2).^2;

    % Calculate the drag force and acceleration due to gravity
    % Drag force (assumes velocity is y(2))
    dydt = [y(2); -DESIGN.g + (drag_force+3*D_airbrake*sin(pi/4)^2) ./ DESIGN.m]; % [dy/dt, dv/dt]
end

function [position,isterminal,direction] = y1_free(t,y)
[T, a, P, rho, nu, mu] = atmosisa(y(1), 'extended', true);

% position = y(2) + 200; % The value that we want to be zero (speed)
    q_current = 0.5 * rho * y(2)^2;
    position = q_current - 10; % Stop when q = 1700 Pa
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

function [t,y,L_D,ic,te,cL] = sadglide(ic,DESIGN,ii)
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
        cL(iter,1) = C_L;
        else
        L_D(iter,1) = L/D; 
        y(iter,1) = y(iter-1,1) + xdistance_traveled;
        y(iter,2) = y(iter-1,2) - ydistance_traveled;
        y(iter,3) = V_Cx;
        y(iter,4) = V_Cy;
        t(iter,1) = iter;
        cL(iter,1) = C_L;
    end
    end

    ic = [y(iter,1), y(iter,2), y(iter,3), y(iter,4)];
    te = t(iter);
end

function [position,isterminal,direction] = y1_glide(t,y)
    position = y(2); % The value that we want to be zero (altitude)
    isterminal = 1;  % Halt integration 
    direction = 0;   % The zero can be approached from either direction
end

function [t,y,L_D,C_L, moment] = funsteady(ic,DESIGN,ii)
% Glide has constant q to keep aero forces
alt = ic(2);

iter = 0;
y(1,1) = ic(1);
y(1,2) = ic(2);
y(1,3) = ic(3);
y(1,4) = ic(4);

while alt > 1700 % && iter < max_iter
    iter = iter + 1; % One iteration is 5 second

     [T, a, P, rho, nu, mu] = atmosisa(alt, 'extended', true);

     if iter == 1
     % finding set q value
     q_con = (0.5*rho*y(iter,3)^2)/2;
      
     end
    
     V = sqrt(2*q_con/(rho));
     W = DESIGN.m*DESIGN.g;
     
    cL = W / (0.5 * rho * V^2 * DESIGN.S(ii));
   

    % Drag Buildup for changing Re and flight condition
    Re = rho.*DESIGN.c.*V./mu;
    [cD0,cDi] = drag(Re,cL,DESIGN,ii);
    cD = (cD0 + cDi);
    
    C_d_airbrake = 0.3;
    D_airbrake = 0.5 * rho * V^2 * DESIGN.S_airbrake*(2/3) * C_d_airbrake;

    L_Dtemp = cL/cD;
    gamma = atan(cD/cL);
    moment(iter,1) = D_airbrake * DESIGN.airbrake_ac;

    if iter == 1
        L_D(iter,1) = L_Dtemp;
        cL(iter,1) = cL;
    else
        y(iter,3) = V*cos(gamma);
        y(iter,4) = V*sin(gamma);
        y(iter,1) = y(iter-1,1) + y(iter,3) * 5;
        y(iter,2) = y(iter-1,2) - y(iter,4) * 5;
        t(iter,1) = iter * 5;
        L_D(iter,1) = L_Dtemp;
        C_L(iter,1) = cL;
    end
    alt = y(iter,2);
end
 
end
