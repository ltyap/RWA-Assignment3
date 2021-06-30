%% steady calculations with flap
clc; clear all; close all;
Uinf = 1;
Winf = 0;
Qinf = sqrt(Uinf^2+Winf^2);
rho = 1.255;
t = 0;
omega = 0;
flap = 1;   % bool property
flap_c = 0.2;     % x/c of flap wrt the chord

styles = {'-', '--', '-.',':', '-','--'};
colors = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.4940, 0.1840, 0.5560];[0.4660, 0.6740, 0.1880]]; % define color scheme, more colors need to be added for larger analysis

% discretization
Npan = 10;
chord = 1;

% for movement to the left at speed Uinf and pitch osciallation at
% frequency omega
kinematics.X0 = -Uinf*t;   % position of LE
kinematics.X0dot = -Uinf;  % X speed
kinematics.Z0 = 0;     % position of LE
kinematics.Z0dot = 0;  % Z speed
kinematics.thetadot = 0;

%% Cl vs alpha
figure()
hold on
it=0;
for beta=deg2rad([0,10,20,30])  % flap angle
    it=it+1;
    Cl = [];
    
    panels = plate_discretization(chord, Npan, flap, beta, flap_c);
    
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
    
    for theta=deg2rad([-15:1:15])
        kinematics.theta = theta;
        temp = [cos(kinematics.theta), sin(kinematics.theta); -sin(kinematics.theta), cos(kinematics.theta)]*[panels.x;panels.z]+[kinematics.X0;kinematics.Z0] ;    %location of TE at t=dt
        X = temp(1,:);
        Z = temp(2,:);
        
        [temp] = LumpedVortex(1, Npan, kinematics, coeff, panels, chord, Qinf, rho);
        Cl = [Cl,temp];
    end
    plot([-15:1:15], Cl, 'linestyle',styles{it}, 'Color', colors(it,:), 'DisplayName',sprintf('beta =%.0f', rad2deg(beta)))
end
grid on
grid minor
legend show
xlabel("\alpha [^\circ]")
ylabel("C_L")

%% Velocity and pressure field
% p_atm = 0;%101300;
% theta = deg2rad(10);
% Ngrid = 100;
%
% [x_grid, z_grid] = meshgrid(linspace(-chord,2*chord,Ngrid),linspace(-1.5*chord,1.5*chord,Ngrid)) ;
% kinematics.theta = theta;
%
% [Cl,dp,~,gamma] = LumpedVortex(1, Npan, kinematics, coeff, panels, chord, Qinf, rho);
% plate_vert.gamma = gamma;   % vertices on the panels
% plate_vert.coord = [cos(kinematics.theta), sin(kinematics.theta); -sin(kinematics.theta), cos(kinematics.theta)]*[panels.xj;panels.zj]+[kinematics.X0;kinematics.Z0] ;
%
% [u_grid, v_grid] = velocity_field(Uinf, Winf, x_grid, z_grid, plate_vert); % calculate velocity in field
% p_grid = pressure_field(Uinf, p_atm, rho, x_grid, z_grid, u_grid, v_grid);
%
% plate = [cos(kinematics.theta), sin(kinematics.theta); -sin(kinematics.theta), cos(kinematics.theta)]*[panels.x;panels.z]+[kinematics.X0;kinematics.Z0] ;
%
% figure()
% hold on
% h = pcolor(x_grid,z_grid,p_grid);
% set(h, 'EdgeColor', 'none');
% colormap pink(32)
% g=colorbar;
% caxis([-1 1])
% shading interp
% % contourf(x_grid,z_grid,p_grid,'Linecolor', 'none')
% plot(plate(1,:),plate(2,:), 'k','LineWidth', 4)
% hold off
% colorbar
% xlim([-chord 2*chord])
% ylim([-1.5*chord 1.5*chord])
% xlabel("x")
% ylabel("z")
% title("Pressure field")
%
% figure()
% hold on
% h = pcolor(x_grid,z_grid,sqrt(u_grid.^2+v_grid.^2));
% set(h, 'EdgeColor', 'none');
% colormap pink(32)
% g=colorbar;
% caxis([0.7 1.5])
% shading interp
% % contourf(x_grid,z_grid,sqrt(u_grid.^2+v_grid.^2),'Linecolor', 'none')
% % colorbar
% t=16;
% quiver(x_grid(1:t:end),z_grid(1:t:end),u_grid(1:t:end),v_grid(1:t:end), 1, 'k')
% plot(plate(1,:),plate(2,:),'k','LineWidth', 4)
% hold off
% xlim([-chord 2*chord])
% ylim([-1.5*chord 1.5*chord])
% xlabel("x")
% ylabel("z")
% title("Velocity field")

