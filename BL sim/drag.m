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
    cDi = cL^2./(DESIGN.cLalpha);

end