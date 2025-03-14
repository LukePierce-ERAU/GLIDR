%% Full simulation of fall from 36.5km for DCR configurations
clear
close all
%% GEOMETRY ASSUMPTION SECTIONS AND GLOBAL VARIABLES ++++++++++++++++++++++


config = 3; % Change to run sims on each config
% 1 ==== Blended Wing Body
% 2 ==== Rogallo Wing
% 3 ==== Bullet Bill


DESIGN = configuration(config);
master = struct();

end_sim = size(DESIGN.S,2);
for ii = 1:1:end_sim

% Freefall until desired q
options = odeset('Events', @y1_free);
place_free = @(t,y) funfree(t,y,DESIGN,ii);
[t,y,te,ye,ei] = ode45(place_free, [0 10000], [36000;0], options); % need to impliment aoa for linear region

for k = 1:numel(t)
    [~,moment(k,:),q(k,:),Re(k,:)] = place_free(t(k),y(k,:));
    % [~,q(k,:)] = place_free(t(k),y(k,:));

end

m_name = ['moment', num2str(ii)];
master.(m_name) = moment;

q_name = ['q', num2str(ii)];
master.(q_name) = q;

Re_name = ['Re', num2str(ii)];
master.(Re_name) = Re;

freevar(:,:,ii) = {t;y(:,1);abs(y(:,2))};
free_out(:,:,ii) = [te,ye,ei];

global pull_init
pull_init = te;

glide_inip(1,ii) = te;

ic = [0;ye(1,1);abs(ye(1,2)); 0];

clear velo

end


%% Cool Plots

free = struct();

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

configs = {'free', 'pull', 'glide', 'steady'};

for jj = 1:1:ii
master.(['t', num2str(jj)]) = vertcat(free.(['t',num2str(jj)]));  %,steady.(['t',num2str(jj)]));
master.(['x', num2str(jj)]) = vertcat(free.(['x',num2str(jj)])); %steady.(['x',num2str(jj)]));
master.(['alt', num2str(jj)]) = vertcat(free.(['alt',num2str(jj)])); %,steady.(['alt',num2str(jj)]));
master.(['speed', num2str(jj)]) = vertcat(free.(['speed',num2str(jj)]));  %steady.(['speed',num2str(jj)]));
end

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
plot(master.moment1, master.alt1, 'b', 'LineWidth', 2);
xlabel('Moment','FontSize', 16);
ylabel('Altitude (m)','FontSize' ,16);
title('Moment vs. Altitude', 'FontSize',18);
grid on;

%% Figure: Airbreak moment vs. Altitude
figure;
plot(master.q1, master.alt1, 'b', 'LineWidth', 2);
xlabel('Dynamic Pressure','FontSize', 16);
ylabel('Altitude (m)','FontSize' ,16);
title('Dynamic Pressure vs. Altitude', 'FontSize',18);
grid on;

%% Figure: Reynolds Number vs. Altitude
figure;
plot(master.Re1, master.alt1, 'b', 'LineWidth', 2);
xlabel('Reynolds Number','FontSize', 16);
ylabel('Altitude (m)','FontSize' ,16);
title('Reynolds Number vs. Altitude', 'FontSize',18);
grid on;

%% Functions for simulations

% Freefall section ode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [dydt,moment,q_save,Re_save] = funfree(t,y,DESIGN,ii)

% Get atmospheric properties at current height
    [T, a, P, rho, nu, mu] = atmosisa(y(1), 'extended', true);

    q = 0.5 * rho * y(2)^2; % Dynamic pressure
    q_save(numel(t),1) = q;

    % Drag Buildup for changing Re and flight condition
    Re = rho.*DESIGN.c.*abs(y(2))./mu;
    Re_save(numel(t),1) = Re;

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
    % q_current = 0.5 * rho * y(2)^2;
    % position = q_current - 10; % Stop when q = 1700 Pa
    target_altitude = 1700;
    position = y(1) - target_altitude ;
    isterminal = 1;  % Halt integration 
    direction = 0;   % The zero can be approached from either direction
end
