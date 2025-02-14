function [Cm0 Cmalpha CmiH CmdE] = cm(DESIGN,cL,)

% % % Blended Wing (BW) Coeff moment (Cm)
% BW.b = 0.45; % wing width
% BW.S = (linspace(1,2,10))'; %area
% BW.AR = BW.b.^2./BW.S; % aspect ratio
% BW.e = 0.35; % efficiency factor
% BW.cLalpha = 1./(pi*BW.e*BW.AR); % lift coefficient relative to AOA
% BW.c = 1.2; % chord length
% BW.M_arm = 0.1; %moment arm length
% 
% 
% BW.Cm0 = 0; % Cm0 = 0 b/c airfoil symmetrical
% 
% BW.Cmalpha = BW.cLalpha.*(-BW.M_arm);
% 
% %% 
% % % % % % % % % % % ROGALLO WING % % % % % % % % % % % % % % % % % % %
% 
% DESIGN.b = 1;
% DESIGN.S = .72; 
% DESIGN.AR = 2; 
% DESIGN.e =.8; 
% DESIGN.cLalpha =1./(pi*DESIGN.e*DESIGN.AR);
% DESIGN.c =1.2 ; 
% DESIGN.M_arm = 0.08;
% DESIGN.M_arm =DESIGN.cLalpha.*(-DESIGN.M_arm);
% 
% 
% % % Horizontal wing
% 
% DESIGN.cLaH = .707;
% DESIGN.nH = 1;
% DESIGN.SH = .6;
% DESIGN.de = 2/3;
% DESIGN.M_armH = .048;
% 
% 
% DESIGN.CmAC = 0 ; % 
% DESIGN.CLO = .3 ; % 
% DESIGN.tauE = .5;


Cm0 = DESIGN.CmAC + (DESIGN.cL).*(-DESIGN.M_arm);

Cmalpha = DESIGN.cLalpha.*(-DESIGN.M_arm)-DESIGN.cLaH.*DESIGN.nH.*(DESIGN.SH/DESIGN.S).*(1-DESIGN.de).*DESIGN.M_armH;

CmiH = -DESIGN.cLaH.*(DESIGN.nH).*(DESIGN.SH/DESIGN.S).*DESIGN.M_armH;

CmdE = -DESIGN.cLaH.*(DESIGN.nH).*(DESIGN.SH/DESIGN.S).*DESIGN.M_armH.*DESIGN.tauE;






