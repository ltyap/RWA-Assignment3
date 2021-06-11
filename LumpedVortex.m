function [p, Cl,new_vortex]=LumpedVortex(coord,c,alpha,N, Qinf, rho, cambered, unsteady,old_vortex,dalpha)
%% Lumped-vortex element
% Input:
% coord - (2,n) table with x,z coordinates of the mean camber line
% c - chord length
% alpha - angle of attack in degrees
% N - number of sub-panels, assuming equal x-length
% Qinf - value of freestream velocity [m/s]
% rho - air density [kg/m3]
% cambered - bool property, 0 for flat plate, 1 for cambered airfoil
% unsteady - bool property, 0 for steady, 1 for unsteady
% old vortex - coordinates and circulation of previous vortex
% dalpha - radial difference between current and previous time step

h = 1/N;     % step size for N uniform grids for x
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

if unsteady
    angle = alpha-0.3*dalpha; % radial position of new vortex
    new_vortex.coord(1) = c/4+3/4*c*cos(angle); % x
    new_vortex.coord(2) = -3/4*c*sin(angle); % z
end

%% Calculate influence matrix coefficients
u = zeros(N, N);
w = zeros(N, N);
coeff = zeros(N, N);
for i=1:N
    for j=1:N
        % assuming singular gamma distribution
        u_ind = v_2D_from_vortex(1, [xj(j), zj(j)], [xi(i), zi(i)]); % not sure about this
        u(i,j) = u_ind(1);
        w(i,j) = u_ind(2);
        coeff(i,j) = dot([u(i,j),w(i,j)],n(i,:));     % influence of element j acting on element i
    end
end

if unsteady
    coeff(:,N+1) = zeros(1,N);
    coeff(N+1,:) = ones(N+1,1);
    for i=1:N
        u_ind = v_2D_from_vortex(1, new_vortex.coord, [xi(i), zi(i)]); %???
        u(i,N+1) = u_ind(1);
        w(i,N+1) = u_ind(2);
        coeff(i,N+1)=dot([u(i,N+1),w(i,N+1)],n(i,:));     % influence of element j acting on element i
    end
end


%% Solving the linear equation
% Right-hand-side vector
rhs=-Qinf*sin(alpha+a);
if unsteady
    rhs=[rhs,old_vortex.gamma]; 
end
gamma=coeff\rhs';   % solve linear equation for vorticity distribution gamma

if unsteady
    new_vortex.gamma = gamma(end);
end

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

% figure()
% hold on
% plot(xi,zi)
% plot(xi_trans,zi_trans)
% hold off
end





