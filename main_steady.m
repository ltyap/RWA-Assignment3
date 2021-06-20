% RWA - assignment 3 - steady cases
clc
clear all
close all

%% Part 1 - steady flow at different AoA
Uinf = 10;  % freestream velocity [m/s]
rho = 1.225;    % air density

aoa = deg2rad([-15:1:15]);   % array of considered angles of attack
chord = 1;
coord(:,1) = linspace(0,chord,100);  % flat plate coordinates (x,z)
coord(:,2) = zeros(1,100);
N = 10; % number of segments

for i=1:length(aoa)
    [~,Cl(i)]=LumpedVortex_different(coord,chord, aoa(i),N, Uinf, rho,0);
end

figure()
plot(rad2deg(aoa),Cl)
grid on
grid minor
xlabel("Angle of attack [^\circ]")
ylabel ("Lift coefficient C_L")

%% trying presure and velocity fields
p_atm = 101300;
aoa = deg2rad(0);
N = 40;

[x_grid, z_grid] = meshgrid(linspace(-chord,2*chord,100),linspace(-1.5*chord,1.5*chord,100)) ;
% p_grid = p_atm*ones(size(x_grid));
% p_grid = rand(size(x_grid));
% u_grid = Uinf*ones(size(x_grid));
% v_grid = zeros(size(x_grid));

[p,Cl,vertices] = LumpedVortex_different(coord, chord, aoa,N, Uinf, rho,0);

[u_grid, v_grid] = velocity_field(Uinf, x_grid,z_grid,vertices); % calculate velocity in field
p_grid = pressure_field(Uinf, p_atm, rho, x_grid, z_grid, u_grid, v_grid);

xi = [0,chord/4,chord];
zi = [0,0,0];
zi = (chord/4-xi)*sin(aoa);
for i=1:length(xi)
    if xi(i)<chord/4
        xi(i) = chord/4-(chord/4-xi(i))*cos(aoa);
    else
        xi(i) = chord/4+(xi(i)-chord/4)*cos(aoa);
    end
end

figure()
hold on
contourf(x_grid,z_grid,p_grid,'Linecolor', 'none')
plot(xi,zi, 'k','LineWidth', 1.5)
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
plot(xi,zi,'k','LineWidth', 1.5)
hold off
xlim([-chord 2*chord])
ylim([-1.5*chord 1.5*chord])
xlabel("x")
ylabel("z")
title("Velocity field")
