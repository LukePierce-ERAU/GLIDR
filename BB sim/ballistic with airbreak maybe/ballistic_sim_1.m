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


config = 3; % Change to run sims on each config
% 1 ==== Blended Wing Body
% 2 ==== Rogallo Wing
% 3 ==== Bullet Bill


DESIGN = configuration(config);
steady = struct();
free = struct();
master = struct();

end_sim = size(DESIGN.S,2);
for ii = 1:1:end_sim


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Freefall until desired q
options = odeset('Events', @y1_free);
place_free = @(t,y) funfree(t,y,DESIGN,ii);
[t,y,te,ye,ei] = ode45(place_free, [0 1000], [36000;0], options); % need to impliment aoa for linear region

for k = 1:numel(t)
    [~,moment(k,:)] = place_free(t(k),y(k,:));
end

    t_name = ['t', num2str(ii)];
    x_name = ['x', num2str(ii)];
    alt_name = ['alt', num2str(ii)];
    speed_name = ['speed', num2str(ii)];
    free.(t_name) = t;
    free.(x_name) = zeros(size(t,1),1);
    free.(alt_name) = y(:,1);
    free.(speed_name) = y(:,2);

freevar(:,:,ii) = {t;y(:,1);abs(y(:,2))};
free_out(:,:,ii) = [te,ye,ei];

global pull_init
pull_init = te;

% glide_inip(1,ii) = te;

if isempty(te)
    glide_inip(1,ii) = 0; % Default to zero if no event is detected
else
    glide_inip(1,ii) = te(1); % Assign the first event time as usual
end

% Pullout @ desired q
% options = odeset('Events', @y1_pull);
% place_pull = @(t,y) funpull(t,y,DESIGN,ii);
% [t,y,te,ye,ei] = ode45(place_pull, [0 5000], [0;ye(1,1);ye(1,2);-pi/2], options); % need to impliment theta and Q

% pull(:,:,ii) = {t;y(:,1);y(:,2);abs(y(:,3));y(:,4)};
% pull_out(:,:,ii) = [te,ye,ei];
% 
% glide_initp(1,ii) = te;

% going from fast to ideal cruise

% options = odeset('Events', @y1_glide);
% chat = @(t,y) funglide(t,y,DESIGN,ii);
% [t,y,te,ye,ei] = ode15s(chat, [0 100], [ye(1,1);ye(1,2);abs(ye(1,3)); 0], options); % need to impliment theta and Q

% ic = [0;ye(1,1);abs(ye(1,2)); 0];
%[t,y,L_D,ic,te] = sadglide(ic,DESIGN,ii);

if isempty(ye)
    ic = [0; y(end,1); abs(y(end,2)); 0]; % Use last calculated state
else
    ic = [0; ye(1,1); abs(ye(1,2)); 0]; % Original behavior
end

% %velo = cell(size(t,1),2);
% velo(:,1) = sqrt(y(:,3).^2+y(:,4).^2);
% %t = zeros(size(velo,1),1);
% glide(:,:,ii) = {t;y(:,1);y(:,2);velo(:,1)};

clear velo

% Stead state cruise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,y,L_D,cL] = funsteady(ic,DESIGN,ii);

velo(:,1) = sqrt(y(:,3).^2+y(:,4).^2);
steadyvar(:,:,ii) = {t;y(:,1);y(:,2);velo(:,1);L_D(:,1);cL(:,1)};
% steady_inip(1,ii) = te; 
if isempty(te)
    steady_inip(1,ii) = NaN; % Assign NaN if no event occurs
else
    steady_inip(1,ii) = te(1); % Take only the first event if multiple exist
end

m_name = ['moment', num2str(ii)];
steady.(m_name) = moment;

    t_name = ['t', num2str(ii)];
    x_name = ['x', num2str(ii)];
    alt_name = ['alt', num2str(ii)];
    speed_name = ['speed', num2str(ii)];
    L_D_name = ['LoveD', num2str(ii)];
    cL_name = ['cL', num2str(ii)];
    steady.(t_name) = t + te;
    steady.(x_name) = y(:,1);
    steady.(alt_name) = y(:,2);
    steady.(speed_name) = velo(:,1);
    steady.(L_D_name) = L_D;
    steady.(cL_name) = cL;






% CL DESIGN POINT EQUATION ++++++++++++++++++++++++++++++++++++++++++++++++

% cL_glide = sqrt(3.*cD0.*(pi.*DESIGN.e.*DESIGN.AR));
% V_trim = sqrt(2.*DESIGN.m.*DESIGN.g./(rho.*DESIGN.S(ii)).*sqrt((1/(pi.*DESIGN.e(ii).*DESIGN.AR(ii)))/(3.*CD0)));
% gamma = atan(CD/CL_glide);
% Glide for min sink^^

% Landing (Bleed V for safe landing)

end


%% Cool Plots







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
master.(['speed', num2str(jj)]) = vertcat(abs(free.(['speed',num2str(jj)])),steady.(['speed',num2str(jj)]));
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
plot(master.cL1(2:end),vertcat(steady.alt1(3:end))/1000,'b','LineWidth',2)
plot(master.cL6(2:end),vertcat(steady.alt6(3:end))/1000,'k','LineWidth',2)
plot(master.cL10(2:end),vertcat(steady.alt10(3:end))/1000,'m','LineWidth',2)
grid on
title('Altitude over Target c_{L}', 'FontSize',18)
xlabel('c_{L}', 'FontSize',16)
ylabel('Altitude [km]', 'FontSize',16)

%% Altitude over Time
figure
hold on
plot(master.t1/60,master.alt1(2:end)/1000,'m','LineWidth',2)
plot(master.t6/60,master.alt6(2:end)/1000,'k','LineWidth',2)
plot(master.t10/60,master.alt10(2:end)/1000,'b','LineWidth',2)
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
xlabel('Dynamic Pressure [Pa]');
ylabel('Altitude [m]');
title('Dynamic Pressure vs. Altitude');
grid on;



% for free_ind = 1:size(free,3)
%     figure(i);
%     free_mat_ind = cell2mat(free(:,:,free_ind));
%     t = free_mat_ind(1:size(free_mat_ind,1)/3);
%     alt = free_mat_ind(size(free_mat_ind,1)/3+1:size(free_mat_ind,1)/3*2);
%     speed = free_mat_ind(size(free_mat_ind,1)/3*2+1:end);
%     hold on;
%     plot(t,alt);
%     title("alt vs time")
%     legend();
%     i = i + 1;
%     figure(i);
%     hold on;
%     plot(speed,alt);
%     title("alt vs speed")
%     legend();
%     i = i - 1;
% end
% % for pull_ind = 1:size(pull,3)
% %     figure(i);
% %     pull_mat_ind = cell2mat(pull(:,:,pull_ind));
% %     t = pull_mat_ind(1:size(pull_mat_ind,1)/5);
% %     x = pull_mat_ind(size(pull_mat_ind,1)/5+1:size(pull_mat_ind,1)/5*2);
% %     alt = pull_mat_ind(size(pull_mat_ind,1)/5*2+1:size(pull_mat_ind,1)/5*3);
% %     speed = pull_mat_ind(size(pull_mat_ind,1)/5*3+1:size(pull_mat_ind,1)/5*4);
% %     hold on;
% %     plot(t+pull_inip(pull_ind),alt);
% %     title("alt vs time")
% %     legend();
% %     i = i + 1;
% %     figure(i);
% %     hold on;
% %     plot(speed,alt);
% %     title("alt vs speed")
% %     legend();
% %     i = i + 1;
% %     figure(i);
% %     hold on;
% %     plot(x,alt);
% %     title("alt vs x-dist")
% %     legend();
% %     i = i - 2;
% % end
% % for glide_ind = 1:size(glide,3)
% %     figure(i);
% %     glide_mat_ind = cell2mat(glide(:,:,glide_ind));
% %     t = glide_mat_ind(1:size(glide_mat_ind,1)/4);
% %     x = glide_mat_ind(size(glide_mat_ind,1)/4+1:size(glide_mat_ind,1)/4*2);
% %     alt = glide_mat_ind(size(glide_mat_ind,1)/4*2+1:size(glide_mat_ind,1)/4*3);
% %     speed = glide_mat_ind(size(glide_mat_ind,1)/4*3+1:end);
% %     hold on;
% %     plot(t+glide_inip(glide_ind),alt);
% %     title("alt vs time")
% %     legend();
% %     i = i + 1;
% %     figure(i);
% %     hold on;
% %     plot(speed,alt);
% %     title("alt vs speed")
% %     legend();
% %     i = i + 1;
% %     figure(i);
% %     hold on;
% %     plot(x,alt);
% %     title("alt vs x-dist")
% %     legend();
% %     i = i - 2;
% % end
% for steady_ind = 1:size(steady,3)
%     figure(i);
%     steady_mat_ind = cell2mat(steady(:,:,steady_ind));
%     t = steady_mat_ind(1:size(steady_mat_ind,1)/4);
%     x = steady_mat_ind(size(steady_mat_ind,1)/4+1:size(steady_mat_ind,1)/4*2);
%     alt = steady_mat_ind(size(steady_mat_ind,1)/4*2+1:size(steady_mat_ind,1)/4*3);
%     speed = steady_mat_ind(size(steady_mat_ind,1)/4*3+1:end);
%     hold on;
%     plot(t+steady_inip(steady_ind),alt);
%     title("alt vs time")
%     legend();
%     i = i + 1;
%     figure(i);
%     hold on;
%     plot(speed,alt);
%     title("alt vs speed")
%     legend();
%     i = i + 1;
%     figure(i);
%     hold on;
%     plot(x,alt);
%     title("alt vs x-dist")
%     legend();
%     i = i - 2;
% end


% Plotting tools
% plot(free(:,1,:),free(:,1,:))



%% Functions for simulations


% Freefall section ode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [dydt,moment] = funfree(t,y,DESIGN,ii)


% Get atmospheric properties at current height
    [T, a, P, rho, nu, mu] = atmosisa(y(1), 'extended', true);

    % Compute Dynamic Pressure (q)
    q = 0.5 * rho * y(2)^2; % Dynamic pressure

    % Drag Buildup for changing Re and flight condition
    Re = rho.*DESIGN.c.*abs(y(2))./mu;
    cL = 0;
    [cD0,cDi] = drag(Re,cL,DESIGN,ii);
    cD = cD0 + cDi;

    % Airbrake Contribution (Assumption: Deployed at 90°)
    C_d_airbrake = 2.2; % Given Cd
    D_airbrake = q *  DESIGN.S_airbrake * C_d_airbrake; % Airbrake drag force
    moment(numel(t),1) = D_airbrake * DESIGN.airbrake_ac; % This is for 1 of 3 airbrakes
    % Calculate the total drag force
    drag_force = 0.5 * cD * DESIGN.S(ii) * rho * y(2)^2 + D_airbrake;

    % Calculate the drag force and acceleration due to gravity
    % Drag force (assumes velocity is y(2))
    dydt = [y(2); -DESIGN.g + drag_force ./ DESIGN.m]; % [dy/dt, dv/dt]
end

function [position,isterminal,direction] = y1_free(t,y)

    
    position = y(1) - 30000; % Event triggers when alt = 30000m
    if position > 0 && y(2) < 5 % Ensure event triggers if velocity is very low
        position = 0; 
    end
    isterminal = 1;  % Halt integration 
    direction = 0;   % The zero can be approached from either direction
end




function [t,y,L_D,C_L,moment] = funsteady(ic,DESIGN,ii)
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
     %q_con = (0.5*rho*y(iter,3)^2)/2;
     q_con = 160;

     end
    
     V = sqrt(2*q_con/(rho));
     
     W = DESIGN.m*DESIGN.g;
     
    cL = W / (0.5 * rho * V^2 * DESIGN.S(ii));
   

    % Drag Buildup for changing Re and flight condition
    Re = rho.*DESIGN.c.*V./mu;
    [cD0,cDi] = drag(Re,cL,DESIGN,ii);
    cD = (cD0 + cDi);

    % Airbrake Contribution
    % C_d_airbrake = 0.3; % Estimate based on "turbulent flow around circular arcs"  Jean-Baptiste R. G. Souppez,1,2 PatrickBot,3 and IgnazioMariaViola
    % D_airbrake = 0.5 * rho * V^2 * DESIGN.S_airbrake*(2/3) * C_d_airbrake;

    % L_Dtemp = cL/cD;
    % gamma = atan(cD/cL);
    % Compute total drag and Lift-to-Drag Ratio (L/D)
    D = 0.5 * rho * V^2 * DESIGN.S(ii) * cD;
    L_D(iter,1) = cL / cD;
    gamma = atan(D / W);
    % moment(iter,1) = D_airbrake * DESIGN.airbrake_ac; % This is for 1 of 3 airbrakes

    % Compute new position
    y(iter+1,1) = y(iter,1) + V * cos(gamma);
    y(iter+1,2) = y(iter,2) - V * sin(gamma);
    y(iter+1,3) = V * cos(gamma);
    y(iter+1,4) = V * sin(gamma);
    alt = y(iter+1,2);

    % Store computed variables
    t(iter,1) = iter;
    C_L(iter,1) = cL;
end

%     if iter == 1
%         L_D(iter,1) = L_Dtemp;
%         cL(iter,1) = cL;
%     else
%         y(iter,3) = V*cos(gamma);
%         y(iter,4) = V*sin(gamma);
%         y(iter,1) = y(iter-1,1) + y(iter,3) * 5;
%         y(iter,2) = y(iter-1,2) - y(iter,4) * 5;
%         t(iter,1) = iter * 5;
%         L_D(iter,1) = L_Dtemp;
%         C_L(iter,1) = cL;
%     end
%     alt = y(iter,2);
% end
% 
% if iter >= max_iter
%         warning('funsteady: Max iterations reached, check logic.');
% 
% end

end


