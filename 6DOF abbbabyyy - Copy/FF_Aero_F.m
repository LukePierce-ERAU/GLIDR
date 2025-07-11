function [F_Aero_B,m_B] = FF_Aero_F(DESIGN,velo_comp,Omega,alpha,beta,delta_ab,rho,mu)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Extraction of Design Variables
ab_S = DESIGN.ab_S;
ab_c = DESIGN.ab_c;
ab_b = DESIGN.ab_b;
ab_t = DESIGN.ab_t;
sb_b = DESIGN.FF_sb_b;
sb_c = DESIGN.FF_sb_c;
sb_Sref = DESIGN.FF_sb_Sref;

f_h = DESIGN.f_h;

% Coefficients to WT test
% Everything listed

velo = sqrt(velo_comp(1)^2 + velo_comp(2)^2 + velo_comp(3)^2);

% aoa correction for aero coefficients


    % Compute Dynamic Pressure (q)
    q = 1/2 * rho * velo^2; % Dynamic pressure

    % Drag Buildup for changing Re and flight condition
    Re = rho*velo/mu .* [f_h  ab_t];

    % Controls stuff and then stuff to find C_delta_ab (cL represents lateral
    % forces because aoa=90 degrees in freefall





% Rotational moment Coefficeint building

c_l_beta = -0.5; % [rad^-1] TO BE TESTED IN WT
c_n_beta = -0.05; % [rad^-1] TO BE TESTED IN WT
c_n_delta_ab = 1; % [rad^-1] TO BE TESTED IN WT
c_m_alpha = 0.5; % [rad^-1] TO BE TESTED IN WT

% Rotational damping Derivative 
c_n_r = -1;
c_m_q = -.33;
c_l_p = -.22;

O_p = Omega(1);
O_q = Omega(2);
O_r = Omega(3);


m_C_l = c_l_beta*beta + sb_b/(2*velo)*c_l_p*O_p; % [ ] Weird, may need coefficeitns added
m_C_n = c_n_delta_ab*delta_ab + c_n_beta*beta + sb_b/(2*velo)*c_n_r*O_r;
m_C_m = c_m_alpha*alpha + sb_c/(2*velo)*c_m_q*O_q;



% Lateral Force Coefficient building ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



[cD0,cDi] = drag(Re,c_n_delta_ab*delta_ab,DESIGN);
C_D = cD0 + cDi;

C_Z = -(C_D(1)+3*C_D(2)); % [ ] Coeficient in the X wind axis % + C_Xq system (only to be used with tail in backwash

C_Ybeta = -1; % [rad^-1] NEED TO TEST IN WT (as a function of Re)
C_Zbeta = -1; % [rad^-1] NEED TO TEST IN WT (THESE ARE ROUGH ESTIMATES)

C_Y = C_Ybeta*beta;

C_X = C_Zbeta*alpha; 

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

F_Aero_B = q*(3*sb_Sref) * [C_X C_Y C_Z];
m_B = q*(3*sb_Sref) * [m_C_l*sb_b m_C_m*sb_c m_C_n*sb_b];
end