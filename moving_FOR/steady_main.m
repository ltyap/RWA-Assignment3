clc; clear all; close all;
Uinf = 10;
Winf = 0;
Qinf = sqrt(Uinf^2+Winf^2);
rho = 1.255;
t = 0;
omega = 0;

% discretization
Npan = 10;
chord = 1;
panels = plate_discretization(chord, Npan);

% for movement to the left at speed Uinf and pitch osciallation at
% frequency omega
kinematics.X0 = -Uinf*t;   % position of LE
kinematics.X0dot = -Uinf;  % X speed
kinematics.Z0 = 0;     % position of LE
kinematics.Z0dot = 0;  % Z speed
kinematics.thetadot = 0;

% self-induced velocity
u = zeros(Npan, Npan);
w = zeros(Npan, Npan);
coeff = zeros(Npan, Npan);
for i=1:Npan
    for j=1:Npan
        % assuming singular gamma distribution
        u_ind = v_2D_from_vortex(1, [panels.xj(j), panels.zj(j)], [panels.xi(i), panels.zi(i)]); % not sure about this
        u(i,j) = u_ind(1);
        w(i,j) = u_ind(2);
        coeff(i,j) = dot([u(i,j),w(i,j)],panels.n(i,:));     % influence of element j acting on element i
    end
end

%% Cl vs alpha
Cl = [];
for theta=deg2rad([-15:1:15])
    kinematics.theta = theta;
%     temp = [cos(kinematics.theta), sin(kinematics.theta); -sin(kinematics.theta), cos(kinematics.theta)]*[panels.x;panels.z]+[kinematics.X0;kinematics.Z0] ;    %location of TE at t=dt
%     X = temp(1,:);
%     Z = temp(2,:);
    
    [temp] = LumpedVortex(1, Npan, kinematics, coeff, panels, chord, Qinf, rho);
    Cl = [Cl,temp];
end

figure()
plot([-15:1:15], Cl)
grid on
grid minor
xlabel("\alpha [^\circ]")
ylabel("C_L")

%% Velocity and pressure field 
p_atm = 101300;
theta = deg2rad(15);
Ngrid = 100;

[x_grid, z_grid] = meshgrid(linspace(-chord,2*chord,Ngrid),linspace(-1.5*chord,1.5*chord,Ngrid)) ;
kinematics.theta = theta;

[Cl,dp,~,gamma] = LumpedVortex(1, Npan, kinematics, coeff, panels, chord, Qinf, rho);
plate_vert.gamma = gamma;   % vertices on the panels
plate_vert.coord = [cos(kinematics.theta), sin(kinematics.theta); -sin(kinematics.theta), cos(kinematics.theta)]*[panels.xj;panels.zj]+[kinematics.X0;kinematics.Z0] ;   

[u_grid, v_grid] = velocity_field(Uinf, Winf, x_grid, z_grid, plate_vert); % calculate velocity in field
p_grid = pressure_field(Uinf, p_atm, rho, x_grid, z_grid, u_grid, v_grid);

plate = [cos(kinematics.theta), sin(kinematics.theta); -sin(kinematics.theta), cos(kinematics.theta)]*[panels.x;panels.z]+[kinematics.X0;kinematics.Z0] ;   

figure()
hold on
contourf(x_grid,z_grid,p_grid,'Linecolor', 'none')
plot(plate(1,:),plate(2,:), 'k','LineWidth', 1.5)
hold off
colorbar
xlim([-chord 2*chord])
ylim([-1.5*chord 1.5*chord])
xlabel("x")
ylabel("z")
title("Pressure field")

figure()
hold on
contourf(x_grid,z_grid,sqrt(u_grid.^2+v_grid.^2),'Linecolor', 'none')
colorbar
quiver(x_grid(1:2:end,1:2:end),z_grid(1:2:end,1:2:end),u_grid(1:2:end,1:2:end),v_grid(1:2:end,1:2:end), 1.5, 'k')
plot(plate(1,:),plate(2,:),'k','LineWidth', 1.5)
hold off
xlim([-chord 2*chord])
ylim([-1.5*chord 1.5*chord])
xlabel("x")
ylabel("z")
title("Velocity field")


