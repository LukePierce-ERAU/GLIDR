function dydt = p_glidr(t,y,DESIGN)
% NOMENCLATURE
% velo_B_E__B  =  Linear velocity of vehicle wrt Earth expressed in the
%                 preffered cordinate body cordinate system. 
%                 DOUBLE UNDERSCORE MEANS EXPRESSED IN CERTAIN FRAME 


% y(1:3) position vectors and are not required in the integration process.
% They are only stored for graphical results.

u = y(4);
v = y(5);
w = y(6);

velo = [u;v;w];

quat0 = y(7);
quat1 = y(8);
quat2 = y(9);
quat3 = y(10);

quat = [quat0; quat1; quat2; quat3];

p = y(11); % B_E__B for all rotations
q = y(12);
r = y(13);

Omega = [p; q; r];

% CONSTANTS

g = 9.81;
% Geometry extraction
I11 = DESIGN.I11;
I13 = DESIGN.I13;
I22 = DESIGN.I22;
I33 = DESIGN.I33;
m = DESIGN.m;

% Integration normalization
quat = quat/norm(quat);

% KINEMATIC EQUATIONS (frame conversion setup and angle tracking) ~~~~~~~~~

T_BL = [quat0^2+quat1^2-quat2^2-quat3^2 2*(quat1*quat2 + quat0*quat3) 2*(quat1*quat3 - quat0*quat2)
        2*(quat1*quat2 - quat0*quat3) quat0^2-quat1^2+quat2^2-quat3^2 2*(quat2*quat3 + quat0*quat1)
        2*(quat1*quat3 + quat0*quat2) 2*(quat2*quat3 - quat0*quat1) quat0^2-quat1^2-quat2^2+quat3^2]; %FILL IN


%%
% These lines below make the assumption that body frame is tied directly to
% wing chord of paraglider wing in glide.

% Actually for paraglider section we can include a mounting angle for the
% paraglider wing to define its placement and not have to rotate the
% body chordinate frame when changing flight phases

% In freefall, these angles are tied to old beta (for new alpha) and old
% zeta (for new beta). Ask Thomas for explanation of this!

    alpha = atan2(w,u);
    beta = asin(v/(sqrt(u^2+v^2+w^2)));

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Atmosphere Subroutine

% Get atmospheric properties at current height
    [T, a, P, rho, nu, mu] = atmosisa(-y(3), 'extended', true);

    clear T a P nu
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% PSI CALC FOR CONTROLS ANALYSIS

psi = atan2(2*(quat1*quat2 + quat0*quat3),quat0^2+quat1^2-quat2^2-quat3^2);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% AERO ANALYSIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% FreeFall Phase (Controller Lives inside Aero block, NEVERMIND it is right here! :) )

if alpha < 2*pi/180 && beta < 2*pi/180
delta_glidr_b = glidr_b_controls(DESIGN,psi);
else
    delta_glidr_b = 0;
end

[F_Aero_B,m_B] = glidr_Aero_F(DESIGN,velo,Omega,alpha,beta,delta_glidr_b,rho,mu);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% EULER's EQUATIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

omega_dot_denom = I11*I33-I13^2;

% l_R not included because it is equal to 0 because we have negligable
% spinning parts

p_dot = 1/omega_dot_denom*(((I22*I33-I33^2-I13^2)*r-I13*(I33+I11-I22)*p)*q+I33*m_B(1)-I13*m_B(3));
q_dot = 1/I22*(((I33-I11)*p)*r+I13*(p^2-r^2)+m_B(2));
r_dot = 1/omega_dot_denom*(((-I11*I22+I11^2+I13^2)*p+I13*(I33+I11-I22)*r)*q+I11*m_B(3)-I13*m_B(1));

Omega_B_E__B_dot = [p_dot; q_dot; r_dot];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% NEWTON'S EQUATIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

s_BE_int = T_BL*velo;

u_dot = r*v-q*w + F_Aero_B(1)/m + T_BL(1,3)*g;
v_dot = p*w-r*u + F_Aero_B(2)/m + T_BL(2,3)*g;
w_dot = q*u-p*v + F_Aero_B(3)/m + T_BL(3,3)*g;

% velo_brc_m = [0 -r q; r 0 -p; -q p 0]; % Velocity body rate components
% velo_B_int = F_Aero_B + T_BL*g_m_L - m*velo_brc_m*velo_B_E__B

velo_B_int = [u_dot; v_dot; w_dot];
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% KINEMATICS EQUATIONS (QUATERNIONS)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

quat_brc_m = [0 -p -q -r; p 0 r -q; q -r 0 p; r q -p 0]; % Quaternion body rate components
quat_dot = 0.5 * quat_brc_m*quat;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


dydt = [s_BE_int;velo_B_int;quat_dot;Omega_B_E__B_dot];
end