

dh = 1000;
ind = round(ic(2)/dh);
y = zeros(ind,4);

for jj  = 1:ind
    [T, a, P, rho, nu, mu] = atmosisa(y(jj,2), 'extended', true);
    if jj == 1
        V = sqrt(ic(3)^2+ic(4)^2);
    else
        V = sqrt(y(jj-1,3)^2+y(jj-1,4)^2);
    end
    W = DESIGN.m.*DESIGN.g;
    gamma = -4.*pi/180;

cD = 0.05;
cL = 0.3; % 2.*W.*cos(gamma)./(rho.*V(jj,1)^2.*DESIGN.S(ii));

L_D(1,jj) = cL/cD;
dx(1,jj) = dh.*L_D(1,jj);

if jj == 1
    y(jj,1) = dx(1,jj);
    y(jj,2) = ic(2);
else
y(jj,1) = y(jj-1,1) + dx(1,jj);
y(jj,2) = ic(2) - jj*dh;
end

t(jj,1) = V.*sqrt(dx(1,jj)^2+dh^2);

if jj == 1
    y(jj,3) = ic(3);
    y(jj,4) = ic(4);
else
y(jj,3) = y(jj-1,3) / dx(1,jj);
y(jj,4) = y(jj-1,4) / dh;
end

% V(jj,1) = 


end