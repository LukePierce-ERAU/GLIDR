function [EXTRACT] = FF_EXTRACT(y, DESIGN)
% Extraction of variables used inside ODE function
%   

% State space redefinition
velo_B = y(:,4:6);
u = velo_B(:,1);
v = velo_B(:,2);
w = velo_B(:,3);

quat0 = y(:,7);
quat1 = y(:,8);
quat2 = y(:,9);
quat3 = y(:,10);

Omega_B = y(:,11:13);
p = Omega_B(:,1);
q = Omega_B(:,2);
r = Omega_B(:,3);

% Euler angle extraction

EXTRACT.psi =   atan2(2.*(quat1.*quat2+quat0.*quat3),(quat0.^2 + quat1.^2 - quat2.^2 - quat3.^2));
EXTRACT.theta = asin(-2.*(quat1.*quat3 - quat0.*quat2));
EXTRACT.phi =   atan2(2*(quat2.*quat3 - quat0.*quat1),(quat0.^2 - quat1.^2 - quat2.^2 + quat3.^2));

% Variable Extraction



EXTRACT.alpha = asin(u./(sqrt(u.^2+v.^2+w.^2)));
EXTRACT.beta = asin(v./(sqrt(u.^2+v.^2+w.^2)));

for ii = 1:size(y,1)

EXTRACT.delta_ab = ab_controls(DESIGN,EXTRACT.psi(ii));

[T, a, P, rho, nu, mu] = atmosisa(-y(3), 'extended', true);
[F_Aero_B,m_B] = FF_Aero_F(DESIGN,velo_B(ii,:),Omega_B(ii,:),EXTRACT.alpha(ii),EXTRACT.beta(ii),EXTRACT.delta_ab,rho,mu);

EXTRACT.F_Aero(ii,:) = F_Aero_B;
EXTRACT.m_B(ii,:) = m_B;
end

end