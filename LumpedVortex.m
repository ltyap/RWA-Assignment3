function [xi, xi_trans, zi_trans, p, Cl,circ, vortex_coord]=LumpedVortex(coord,c,alpha,N, Qinf, rho, cambered, unsteady, circ_old, vortex_coord_old, dalpha)
%% Lumped-vortex element
% Input:
% coord - (2,n) table with x,z coordinates of the mean camber line
% c - chord length 
% alpha - angle of attack in degrees
% N - number of sub-panels, assuming equal x-length
% Qinf - value of freestream velocity [m/s]
% rho - air density [kg/m3] 
% cambered - bool property, 0 for flat plate, 1 for cambered airfoil
% unsteady - bool property, 0 for steady calculations, 1 for unsteady
% circ_old - value of circulation for previous time step
% dalpha - radial pitch for one time step

h = 1/N;     % step size for N uniform grids for x  
alpha = alpha/180*pi ;  % translating the angle of attack into radians

%% Discretize camber into points(panel edges)
x=0:h:c;
if cambered
    z = interp1(coord(1,:),coord(2,:),x);
	z(1)=0;
else % flat plate
    z = zeros(1,N+1);
end
xj = zeros(1, N);
xi = zeros(1, N);
zj = zeros(1, N);
zi = zeros(1, N);
a = zeros(1, N);
n = zeros(N, 2);
for i =1:N      % loop over each panel
   % Vortex points
    xj(i)= x(i)+1/4*(x(i+1)-x(i));
    zj(i)= z(i)+1/4*(z(i+1)-z(i));
    
    % Collocation points
    xi(i)= x(i)+3/4*(x(i+1)-x(i)); 
    zi(i)= z(i)+3/4*(z(i+1)-z(i));
    
    % Normal vectors and panel slopes in collocation points
    a(i)=atan2((z(i+1)-z(i)),(x(i+1)-x(i)));    % angle of panel used for calculating normal vector
    n(i,:) = [cos(a(i)), sin(a(i))];
    n(i,:)= ([cos(pi/2), -sin(pi/2); sin(pi/2), cos(pi/2)]*n(i,:)')';   %normal vector of i-th panel
    a(i)=asin(n(i,1));  % angle of panel consistent with the following equations
end

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
        Gamma=1;    % assuming singular gamma distribution
        u_ind = v_2D_from_vortex(Gamma, [xj(j), zj(j)], [xi(i), zi(i)]); % not sure about this
        u(i,j) = u_ind(1);
        w(i,j) = u_ind(2);
%         r_sq=(xi(i)-xj(j))^2 + (zi(i)-zj(j))^2;
%         u(i,j)=Gamma/(2*pi*r_sq)*(zi(i)-zj(j));     % velocity of arbitrary point (x-coordinate)
%         w(i,j)=-Gamma/(2*pi*r_sq)*(xi(i)-xj(j));    % velocity of arbitrary point (z-coordinate)
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
rhs=-Qinf*sin(alpha+a);
if unsteady
    rhs=[rhs,circ_old]; 
end
gamma=coeff\rhs';   % solve linear equation for vorticity distribution gamma

%% Computing pressure and lift
p = zeros(1, length(gamma));
dL = zeros(1, length(gamma));
for i=1:N   % loop over each panel
    p(i)=rho*Qinf*gamma(i)/sqrt((x(i+1)-x(i))^2+(z(i+1)-z(i))^2);  % pressure difference on panel
    dL(i)=rho*Qinf*gamma(i);    % lift difference of panel
end

dcp = p/(0.5*rho*Qinf^2*c);     % pressure coefficient distribution
L=sum(dL);      % total lift
Cl=L/(0.5*rho*Qinf^2*c);    % total lift coefficient

x0 = c/4;% quarter chord points around which we rotate flat plate
z0 = 0;
xi_trans = zeros(1,N);
zi_trans = zeros(1,N);
if ~cambered % translate coordinates of collocation points
    xi_trans = xi*cos(alpha);
    zi_trans = (x0-xi)*sin(alpha);
end

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





