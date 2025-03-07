function [cD0,cDi] = drag(Re,cL,DESIGN,ii)
%UNTITLED3 Drag buildup at a given speed and altitude
%   based on atmospheric conditions a drag buildup is made to estimate drag
%   for all designs for GLIDR Capstone. Method is based on general skin
%   friction coefficent estimates. To be expanded on as model becomes more
%   sophisticated.
   
% To be removed with global geometry variables

   
if Re < 1000
    Re = 1000;
end
   
    if Re < 5e5
        cf = 1.328/sqrt(Re);
    else
        cf = 0.074/Re^0.2;
    end
    cD0 = 1.25.*cf.*DESIGN.S_ratio;
    % cDi = cL^2./(DESIGN.cLalpha(ii));
    
    if length(DESIGN.cLalpha) > 1
    cDi = cL^2 ./ DESIGN.cLalpha; % ✅ Works if `DESIGN.cLalpha` is a vector
else
    cDi = cL^2 / DESIGN.cLalpha; % ✅ Works if `DESIGN.cLalpha` is a scalar
    end

% % Airbrake drag calculation
% A_brake = 0.05;  % Airbrake frontal area in m² (example value, adjust as needed)
% Cd_brake = 2.2;  % Given drag coefficient for worst-case scenario
% Fd_airbrake = q * A_brake * Cd_brake;


end