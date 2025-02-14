% Aero analysis Assumptions for BL vehicle
    S_ratio = 2.5; % Made from similar aircraft; original eqation: S_wet/S_ref
    AR = 80/120; % Hard to guess right now
    e = 0.35; % Low Estimate after looking at sources
    cLalpha = 1/(pi*e*AR);
    S = 1.75; % Very crude guess based on prelim CAD work




% ODE for position and velocity
options = odeset('Events', @y1_0);
[t,y,te,ye,ei] = ode45(@(t,y) vap1(t,y), linspace(0,200,1000), [36000;0;-pi/2;0], options);

function dydt = vap1(t, y)
    g = 9.81; % ft/sec^2
    m = 8; % Mass in kg 
    c = 1.2; % meters

    % Get atmospheric properties at current height
    [T, a, P, rho, nu, mu] = atmosisa(y(1), 'extended', true);


% Drag Buildup for changing Re and flight condition
    
    Re = rho.*c.*y(2)./mu;
if Re < 1000
    Re = 1000;
end
    

    if y(1)>30000
        cL = 0;
    else
        cL = 0.5;
    end


    
    if Re < 5e5
        cf = 1.328/sqrt(Re);
    else
        cf = 0.074/Re^0.2;
    end
    cD0 = 1.25*cf*S_ratio;
    
    cD = cD0 + cL^2/(cLalpha);
    LD = cL/cD;


   
    cD = cD0 + cL^2/cLalpha;

    
    drag_force = 0.5 * cD * S * rho * y(2)^2;


    % Calculate the drag force and acceleration due to gravity
     % Drag force (assumes velocity is y(2))
    dydt = [y(2); -g + drag_force / m]; % [dy/dt, dv/dt]
end

function [position,isterminal,direction] = y1_0(t,y)
    position = y(1); % The value that we want to be zero (altitude)
    isterminal = 1;  % Halt integration 
    direction = 0;   % The zero can be approached from either direction
end


[T_f, a_f, P_f, rho_f,nu_f,mu_f] = atmosisa(y(:,1),extended=true);

Mach = abs(y(:,2))./a_f;

figure
plot(Mach,y(:,1))
xlabel('Mach number')
ylabel('Altitude in ft')
grid on

