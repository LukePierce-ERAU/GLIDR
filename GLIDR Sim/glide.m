function [V_opt,dydt,dx,dt] = glide(cL,cD,DESIGN,ii)




dydt = sqrt(2*DESIGN.m.*DESIGN.g./(rho.*DESIGN.S(ii))).*cD./cL^(3/2);



end