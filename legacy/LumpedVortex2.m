function [Cl,p]=LumpedVortex2(xi,zi,xj,zj,n,c,unsteady,Qinf,rho,alpha)
%[xi, xi_trans, zi_trans, p, Cl,circ, vortex_coord]=LumpedVortex(coord,c,alpha,N, Qinf, rho, unsteady, circ_old, vortex_coord_old, dalpha)
%% Lumped-vortex element
% Input:
% coord - (2,n) table with x,z coordinates of the mean camber line
% c - chord length
% alpha - angle of attack in radians
% N - number of sub-panels, assuming equal x-length
% Qinf - value of freestream velocity [m/s]
% rho - air density [kg/m3]
% unsteady - bool property, 0 for steady calculations, 1 for unsteady
% circ_old - value of circulation for previous time step
% dalpha - radial pitch for one time step
N = length(xi);

%% Calculate influence matrix coefficients
if unsteady     % additional DOF for the unknown wake
    u = zeros(N+1, N+1);
    w = zeros(N+1, N+1);
    coeff = zeros(N+1, N+1);
else % steady calculations
    u = zeros(N, N);
    w = zeros(N, N);
    coeff = zeros(N, N);
end

for i=1:N
    for j=1:N
        Gamma = 1;    % assuming singular gamma distribution
        u_ind = v_2D_from_vortex(Gamma, [xj(j), zj(j)], [xi(i), zi(i)]); % not sure about this
        u(i,j) = u_ind(1);
        w(i,j) = u_ind(2);
        coeff(i,j)=dot([u(i,j),w(i,j)],n(i,:));     % influence of element j acting on element i
    end
end

if unsteady
    coeff(N+1,:) = ones(N+1,1);
    for i=1:N
        u_ind = v_2D_from_vortex(Gamma, vortex_coord_old, [xi(i), zi(i)]); %???
        u(i,N+1) = u_ind(1);
        w(i,N+1) = u_ind(2);
        coeff(i,N+1)=dot([u(i,N+1),w(i,N+1)],n(i,:));     % influence of element j acting on element i
        %         coeff(:,i) = ; % influence of wake vortex
    end
end

%% Solving the linear equation
% Right-hand-side vector
rhs=-Qinf*sin(alpha)*ones(N,1);
if unsteady
    rhs=[rhs;circ_old];
end
gamma=coeff\rhs;   % solve linear equation for vorticity distribution gamma

%% Computing pressure and lift
p = zeros(1, length(gamma));
dL = zeros(1, length(gamma));
for i=1:N   % loop over each panel
    %     p(i)=rho*Qinf*gamma(i)/sqrt((x(i+1)-x(i))^2+(z(i+1)-z(i))^2);  % pressure difference on panel
    dL(i)=rho*Qinf*gamma(i);    % lift difference of panel
end

% dcp = p/(0.5*rho*Qinf^2*c);     % pressure coefficient distribution
L=sum(dL);      % total lift
Cl=L/(0.5*rho*Qinf^2*c);    % total lift coefficient

% x0 = c/4;% quarter chord points around which we rotate flat plate
% z0 = 0;
xi_trans = zeros(1,N);
zi_trans = zeros(1,N);
% if ~cambered % translate coordinates of collocation points
%     xi_trans = x0+xi*cos(alpha);
%     zi_trans = -xi*sin(alpha);
% end

% figure()
% hold on
% plot(xi,zi)
% plot(xi_trans,zi_trans)
% hold off
%% Calculate circulation and its coordinates
if unsteady
    circ = gamma(end);
    %     vortex_coord = [1,1]; along the path of the trailing edge, 0.2-0.3 of
    %     the movement in the latest time step
    %   only for pitching movement
    vortex_coord(1) = 1/4*c+3/4*c*cos(alpha-dalpha);
    vortex_coord(2) = -3/4*c*sin(alpha-dalpha);
end

end