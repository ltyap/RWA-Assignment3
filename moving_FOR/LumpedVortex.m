function [Cl,dp, vertices,gamma, Cl_c, Cl_nc] = LumpedVortex(steady, N, kinematics, coeff, panels, chord, Qinf, rho, vertices, TE_old, TE_new, gamma_old, dt)
temp = [cos(kinematics.theta), sin(kinematics.theta); -sin(kinematics.theta), cos(kinematics.theta)]*...
    [panels.xi;panels.zi]+[kinematics.X0; kinematics.Z0];
panels.Xi = temp(1,:);
panels.Zi = temp(2,:);

if ~steady  % coordinates of new vortex - for now just the avergae TO DO!!!!
    new_vortex.x = 0.3*TE_old(1)+0.7*TE_new(1);
    new_vortex.z = 0.3*TE_old(2)+0.7*TE_new(2);
end

temp = [cos(kinematics.theta), -sin(kinematics.theta); sin(kinematics.theta), cos(kinematics.theta)]*...
    [-kinematics.X0dot;-kinematics.Z0dot]+[-kinematics.thetadot*panels.eta;kinematics.thetadot*panels.xi-panels.deta];
U = temp(1,:);%*ones(1,N);
W = temp(2,:);%*ones(1,N);

if ~steady
    % velocities induced by previous wakes on control points
    for i=1:length(panels.xi)
        u_temp = 0;
        w_temp = 0;
        for j=1:length(vertices.gamma)
            u_ind = v_2D_from_vortex(vertices.gamma(j), [vertices.x(j), vertices.z(j)], [panels.Xi(i), panels.Zi(i)]);
            u_temp = u_temp+u_ind(1);
            w_temp = w_temp+u_ind(2);
        end
        u_w(i) = u_temp;
        w_w(i) = w_temp;
        % velocity induced by unknown circulation
        u_ind = v_2D_from_vortex(1, [new_vortex.x, new_vortex.z], [panels.Xi(i), panels.Zi(i)]);
        coeff(i,N+1)=dot(u_ind,panels.n(i,:));
    end
    coeff(N+1,:) = ones(1,N+1);
end

% right-hand side
if steady
    rhs = -dot([U; W],panels.n');
else
    rhs = -dot([U+u_w; W+w_w],panels.n');
    rhs = [rhs, sum(gamma_old)];
end

% solve system
gamma = coeff\rhs';

%change!!!
% panels.a = kinematics.theta*ones(N,1);  % angle of panel consistent with the following equations
% panels.tau = [cos(panels.a), -sin(panels.a)]; % tangential components 

% loads
dp = zeros(1, length(gamma));
dL = zeros(1, length(gamma));
for i=1:N   % loop over each panel
    dl(i) = sqrt((panels.x(i)-panels.x(i+1))^2+(panels.z(i)-panels.z(i+1))^2);    % panel length
    if steady
        dp(i) = rho*(dot([U(i);W(i)],panels.tau(i,:))*gamma(i)/dl(i));   % pressure difference on panel
    else
        dgamma = sum(gamma(1:i))-sum(gamma_old(1:i));
        dp(i) = rho*(dot([U(i)+u_w(i);W(i)+w_w(i)],panels.tau(i,:))*gamma(i)/dl(i)+dgamma/dt);   % pressure difference on panel
    end
    dL(i) = dp(i)*dl(i)*cos(panels.a(i));    % lift difference of panel
    dL_c(i) = rho*Qinf*gamma(i);    % circulatory
end

dcp = dp/(0.5*rho*Qinf^2*chord);     % pressure coefficient distribution
L = sum(dL);      % total lift - circulatory+uncirculatory
L_c = sum(dL_c);    % total lift - only circulatory (so from Kutta condition)
L_nc = L-L_c;   % lift - non-circulatory
Cl_c  = L_c/(0.5*rho*Qinf^2*chord); % circulatory loads
Cl_nc  = L_nc/(0.5*rho*Qinf^2*chord); % non-circulatory loads
Cl = L/(0.5*rho*Qinf^2*chord);    % total lift coefficient

% add new vortex
if ~steady
    vertices.x = [vertices.x,new_vortex.x];
    vertices.z = [vertices.z,new_vortex.z];
    vertices.gamma = [vertices.gamma,gamma(end)];
else
    vertices = [];
end

% Qt = dot([U+u_w',W+w_w'],panels.tau)

% dL=rho*Qinf*gamma; % bad hack
% CL = sum(dL)/(0.5*rho*Qinf^2*chord);

end