% 2D velocity induced by a vortex filament
function u_ind = v_2D_from_vortex(gamma, vp, cp)
% point1 = point1(:);
% point2 = point2(:);
% cp = cp(:);
r = sqrt(sum((cp-vp).^2));

XVP = vp(1);
ZVP = vp(2);

XCP = cp(1);
ZCP = cp(2);
% induced velocity
u_ind = gamma/(2*pi*r^2)*[0,1;-1,0]*[XCP-XVP;ZCP-ZVP];
end