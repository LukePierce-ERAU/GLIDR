% Main ode45 6-DOF GLIDR simulation
clear
close all

% Design definition function
[DESIGN] = BP_geo();



% Simulation ic

s_ic_B = [0; 0; -36500]; % m
velo_B_ic = [0.2; -01; 5];
quat_ic = [1; 0; 0; 0];
Omega_ic = [0.5; -0.4; 0];


%% Simulation Execution
% Multiple simulations are excecuted because of the harsh change in
% aerodynamic environment. Breaking it up like helps helps ode45 and allows
% for a non continuios (<--- I tried on the spelling) derivative (less hard for computer to estimate)

% Freefall Phase
e_phase = @(t, y) phase(t, y, DESIGN);
options1 = odeset('Events', e_phase);

% Debug tracking process


[t1,y1,FORCES,te,ye,ie,] = ode15s(@(t,y)free_fallin(t,y,DESIGN),[0 15],[s_ic_B;velo_B_ic;quat_ic;Omega_ic],options1);

[EXTRACT] = FF_EXTRACT(y1,DESIGN);

quat = y1(:,7:10);

% Paraglide Phase
% e_term = @(t, y) ground(t, y, DESIGN);
% options2 = odeset('Events', e_term);
% 
% [t2,y2] = ode45(@(t,y)p_glidr(t,y,DESIGN),[0 50],ye(1:3),ye(4:6),ye(7:10),ye(11:13),options2);
% 
% t = [t1;t2];
% y = [y1;y2];

% Data Processing including conversion from quaternions to Euler angles.
% And then graphs :)

% And then bug jason for animation....

figure
plot(y1(:,1),y1(:,2))
grid on
title('Postition')

figure
plot(t1,-y1(:,3))
grid on
title('Alt vs time')

figure
hold on
plot(t1,y1(:,4))
plot(t1,y1(:,5))
plot(t1,y1(:,6))
grid on
title('velocity vs time')
legend('u','v','w')

figure
hold on
plot(t1,y1(:,11))
plot(t1,y1(:,12))
plot(t1,y1(:,13))
title('Omega vs time')
legend('p','q','r')
grid on

figure
hold on
plot(t1,EXTRACT.theta./(2*pi))
plot(t1,EXTRACT.phi./(2*pi))
plot(t1,EXTRACT.psi./(2*pi))
title('Euler angles vs time')
legend('theta','phi','psi')
grid on

figure
hold on
plot(t1,EXTRACT.F_Aero(:,1))
plot(t1,EXTRACT.F_Aero(:,2))
plot(t1,EXTRACT.F_Aero(:,3))
title('Aero forces vs time')
legend('x','y','z')
grid on

figure
hold on
plot(t1,EXTRACT.alpha)
plot(t1,EXTRACT.beta)
title('Aero angles vs time')
legend('alpha','beta')
grid on