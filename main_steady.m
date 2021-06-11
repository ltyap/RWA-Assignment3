% RWA - assignment 3 - steady cases
clc
clear all
close all

%% Part 1 - steady flow at different AoA
Uinf = 10;  % freestream velocity [m/s]
rho = 1.225;    % air density

aoa = [-15:1:15];   % array of considered angles of attack
chord = 1;
coord(:,1) = linspace(0,chord,100);  % flat plate coordinates (x,z)
coord(:,2) = zeros(1,100);
N = 10; % number of segments

for i=1:length(aoa)
    [~,~,~,~,Cl(i)]=LumpedVortex(coord,chord, aoa(i),N, Uinf, rho,0);
end

figure()
plot(aoa,Cl)
grid on
grid minor
xlabel("Angle of attack [^\circ]")
ylabel ("Lift coefficient C_L")

%% trying presure and velocity fields
p_atm = 101300;
aoa=15;
N=10;

[x_grid, z_grid] = meshgrid(linspace(-chord,2*chord,3*N),linspace(-1.5*chord,1.5*chord,3*N)) ;
% p_grid = p_atm*ones(size(x_grid));
p_grid = rand(size(x_grid));
u_grid = Uinf*ones(size(x_grid));
v_grid = zeros(size(x_grid));

[~,xi, zi,p,Cl]=LumpedVortex(coord,chord, aoa,N, Uinf, rho,0);

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
quiver(x_grid,z_grid,u_grid,v_grid, 0.5)
plot(xi,zi,'k','LineWidth', 1.5)
hold off
xlim([-chord 2*chord])
ylim([-1.5*chord 1.5*chord])
xlabel("x")
ylabel("z")
title("Velocity field")
