%% Main
clear
close all
clc
%% Read data
naca0015 = fscanf(fopen('Cp_NACA0015.txt', 'r'), '%g %g %g %g', [4 160]);
naca0015_a0_up = flip(naca0015(1:2, 1:80), 2);
naca0015_a0_low = naca0015(1:2, 81:end);
naca0015_a5_up = flip(naca0015(3:4, 1:80), 2);
naca0015_a5_low = naca0015(3:4, 81:end);

naca0015_a0 = [flip(naca0015_a0_up, 2), naca0015_a0_low];
naca0015_a5 = [flip(naca0015_a5_up, 2), naca0015_a5_low];
dcp_naca0015_a5 = naca0015_a5_low(2, :) - naca0015_a5_up(2, :);
dcp_naca0015_a0 = naca0015_a0_low(2, :) - naca0015_a0_up(2, :);

coord=[naca0015_a5_up(1, :); zeros(1, length(naca0015_a5_up))];

%% NACA 4-series airfoils coordinates
x = linspace(0, 1, 100); % x coordinates for NACA 4-series airfoils used in camberline.m
z_2412 = camberline(x, 2, 4); % z coordinates for NACA 4-series airfoils

%% Comparison of various methods
N = 100; % number of panels for calculations with LVM 

% NACA 2412
[x_2412_num15, dcp_2412_num15, Cl_2412_num15]=LumpedVortex([x; z_2412],5,15, 1, 1.225); % pressure difference distribution and Cl calculation for NACA2412 using LVM
% [Cp_2412_thin_airfoil_mid, ~, Cl_2412_thin] = ThinAirfoilMidpoints(5, [x; z_2412]); % pressure difference distribution and Cl calculation for NACA2412 using Thin Airfoil Theory

figure
h1 = axes;
hold on
% plot(x, Cp_2412_thin_airfoil_mid)
% plot(x_2412_num5, dcp_2412_num5, 'gx-')
plot(x_2412_num15, dcp_2412_num15, 'x-')
% plot(x_2412_num30, dcp_2412_num30, 'bx-')
% plot(x_2412_num100, dcp_2412_num100, 'rx-')
% plot(x, Cp_2412_thin_airfoil_mid)
% legend('N = 5 {LVM}', 'N = 15 {LVM}', 'N = 30 {LVM}', 'N = 100 {LVM}')
xlabel('x/c')
ylabel('dCp')
ylim([0, 5])
grid on
grid minor
hold off

%% NACA 0015
[x_0015_numm5, dcp_0015_numm5, Cl_0015_numm5]=LumpedVortex(coord,-5,N, 1, 1.225); % pressure difference distribution and Cl calculation for NACA0015 using LVM
% [Cp_0015_thin_airfoil_mid, ~, Cl_0015_thin] = ThinAirfoilMidpoints(5, [x; zeros(1, length(x))]); % pressure difference distribution and Cl calculation for NACA2412 using Thin Airfoil Theory

figure 
h2 = axes;
hold on
% plot(x, Cp_0015_thin_airfoil_mid)
% plot(x_0015_num0, dcp_0015_num0)
% plot(naca0015_a0_up(1, :), dcp_naca0015_a0)
% plot(x_0015_num5, dcp_0015_num5, 'b-x')
% plot(x_0015_num10, dcp_0015_num10, 'g-')
% plot(x_0015_num30, dcp_0015_num30, 'g-x')
% plot(x_0015_num100, dcp_0015_num100, 'm-x')
% plot(naca0015_a5_up(1, :), dcp_naca0015_a5, 'r-', 'LineWidth',2)
% legend('N=5 (LVM)', 'N=30 (LVM)', 'N=100 (LVM)', 'data')
plot((-15:5:15), [Cl_0015_numm15, Cl_0015_numm10, Cl_0015_numm5, Cl_0015_num0, Cl_0015_num5, Cl_0015_num10, Cl_0015_num15])
xlabel('Angle of attack [dgr]')
ylabel('Lift coefficient')
line([0 0], [-2 2], 'Color', 'k');  %x-axis
line([-15 15], [0 0], 'Color', 'k');  %y-axis
xlim([-15, 15])
grid on
grid minor
hold off

%% Example 1 p.268 in "Low-Speed Aerodynamics" by J. Katz and A. Plotkin
alpha = 0;
eps = 0.1;
n = 100;
x = linspace(0, 1, n+1);
z_new = parabCamber(x, eps);
coord = [x; z_new];
[x_ex, dcp_ex, cl_ex] = LumpedVortex(coord,alpha,n);
dCp = zeros(1, length(x));
for i = 1:n
    dCp(i) = 4*sqrt((1-x(i))/x(i))*alpha*pi/180 + 32*eps/1*sqrt(x(i)/1*(1-x(i)/1));
end
figure
hold on 
plot(x, dCp)
plot(x_ex, dcp_ex, 'o')
hold off