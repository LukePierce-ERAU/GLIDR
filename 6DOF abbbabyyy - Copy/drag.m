function [cD0,cDi] = drag(Re,cL,DESIGN)
%UNTITLED3 Drag buildup at a given speed and altitude
%   based on atmospheric conditions a drag buildup is made to estimate drag
%   for all designs for GLIDR Capstone. Method is based on general skin
%   friction coefficent estimates. To be expanded on as model becomes more
%   sophisticated.
   
% To be removed with global geometry variables

for ii = 1:length(Re)
    if Re(ii) < 1000
        Re(ii) = 1000;
    end
   
    if Re(ii) < 5e5
        cf(ii) = 1.328/sqrt(Re(ii));
    else
        cf(ii) = 0.074/Re(ii)^0.2;
    end
end
    cD0 = 1.25*cf.*[DESIGN.ab_Swet  DESIGN.f_Swet];
    % cDi = cL^2./(DESIGN.cLalpha(ii));
    

    cDi = cL.^2 ./ DESIGN.ab_cLalpha; % ✅ Works if `DESIGN.cLalpha` is a scalar
  

% % Airbrake drag calculation
% A_brake = 0.05;  % Airbrake frontal area in m² (example value, adjust as needed)
% Cd_brake = 2.2;  % Given drag coefficient for worst-case scenario
% Fd_airbrake = q * A_brake * Cd_brake;


end