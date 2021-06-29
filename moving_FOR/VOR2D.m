function [u_ind, w_ind] = VOR2D(X,x_vp,Z,z_vp,i,j,gamma)

r = sqrt((X(i)-x_vp(j)).^2 + (Z(i)-z_vp(j)).^2);

XVP = x_vp(j);
ZVP = z_vp(j);

XCP = X(i);
ZCP = Z(i);
% induced velocity
temp = gamma/(2*pi*r^2)*[0,1;-1,0]*[XCP-XVP;ZCP-ZVP];
u_ind = temp(1);
w_ind = temp(2);
end