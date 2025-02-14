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

% Glide @ V for L/D max
options = odeset('Events', @y1_glide);
place_glide = @(t,y) funglide(t,y,DESIGN,ii);
[t,y,te,ye,ei] = ode45(place_glide, [0 100], [ye(1,1);ye(1,2);ye(1,3); 0], options); % need to impliment theta and Q



glide(:,:,ii) = {t;y(:,1);y(:,2);abs(y(:,3))};





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
cD = cD0 + cDi;

drag_force = 0.5 * cD * DESIGN.S(ii) * rho * y(2).^2;

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
        cL = 0.5;
    end
   
    
[cD0,cDi] = drag(Re,cL,DESIGN,ii);
cD = cD0 + cDi;
    
drag_force = 0.5 .* rho .* y(3)^2 .* DESIGN.S(ii) .* DESIGN.S_ratio .* cD;
lift_force = 0.5 .* rho .* y(3)^2 .* DESIGN.S(ii) .* cL;
   
    
% Q assuming aircraft is point mass. Can be updated with inertia. Replace m
% with Izz to do so. ++++++++++++++++++++++++++++++++++++++++++++++++++++

% Assumes that aero center is on same waterline as cg.++++++++++++++++++++

Q = ( .5*rho*(y(3)^2)*cL*(DESIGN.S(ii)/DESIGN.m) - DESIGN.g*cos(y(4)) )/abs(y(3));
    
 % Calculate the drag force and acceleration due to gravity
     % Drag force (assumes velocity is y(2))
    dydt = [abs(y(3)).*cos(y(4)); abs(y(3)).*sin(y(4)); -DESIGN.g.*sin(y(4))+drag_force/DESIGN.m; Q]; % [dx/dt, dy/dt, dv/dt, dtheta/dt]
end

function [position,isterminal,direction] = y1_pull(t,y)
    position = y(4); % The value that we want to be zero (theta)
    isterminal = 1;  % Halt integration 
    direction = 0;   % The zero can be approached from either direction
end

function dydt = funglide(t,y,DESIGN,ii)

% Get atmospheric properties at current height
    [T, a, P, rho, nu, mu] = atmosisa(y(1), 'extended', true);

    % Drag Buildup for changing Re and flight condition
    Re = rho.*DESIGN.c.*y(2)./mu;


    W = DESIGN.m.*DESIGN.g;

cL = 2.*W./(rho.*DESIGN.S(ii).*y(3)^2);
[cD0W,cDiW] = drag(Re,cL,DESIGN,ii);


% V_trim_g = sqrt(2.*DESIGN.m.*DESIGN.g./(rho.*DESIGN.S(ii)).*sqrt((1/(pi.*DESIGN.e(ii).*DESIGN.AR(ii)))/(3.*cD0)));
% cLW = sqrt(3.*cD0.*(pi.*DESIGN.e.*DESIGN.AR));


cLH = DESIGN.S(ii).*cL.*DESIGN.M_arm/(DESIGN.SH.*DESIGN.M_armH);
[cD0H,cDiH] = drag(Re,cLH,DESIGN,ii);

cD = cD0W + cDiW + cD0H + cDiH;
drag_force = 0.5 .* rho .* y(3)^2 .* DESIGN.S(ii) .* DESIGN.S_ratio .* cD;
lift_force = 0.5 .* rho .* y(3)^2 .* DESIGN.S(ii) .* cL;
%V_trim_g = sqrt(2.*DESIGN.m.*DESIGN.g./(rho.*DESIGN.S(ii)).*sqrt((1/(pi.*DESIGN.e(ii).*DESIGN.AR(ii)))/(3.*cD0)));

gamma = abs(atan(y(2)/y(1)))

x_accel = -drag_force.*cos(gamma)/DESIGN.m - lift_force.*sin(gamma)/DESIGN.m;
y_accel =  -drag_force.*sin(gamma)/DESIGN.m + lift_force.*cos(gamma)/DESIGN.m - DESIGN.g;

dydt = [y(3); y(4); -drag_force.*cos(gamma)/DESIGN.m - lift_force.*sin(gamma)/DESIGN.m; -drag_force.*sin(gamma)/DESIGN.m + lift_force.*cos(gamma)/DESIGN.m - DESIGN.g];
% Glide for min sink^^

end


function [position,isterminal,direction] = y1_glide(t,y)
    position = y(2); % The value that we want to be zero (altitude)
    isterminal = 1;  % Halt integration 
    direction = 0;   % The zero can be approached from either direction
end