clear 
close all
%% Loading NACA0015 data points
naca0015 = fscanf(fopen('Cp_NACA0015.txt', 'r'), '%g %g %g %g', [4 160]);
naca0015_a0_up = flip(naca0015(1:2, 1:80), 2);
naca0015_a0_low = naca0015(1:2, 81:end);
naca0015_a5_up = flip(naca0015(3:4, 1:80), 2);
naca0015_a5_low = naca0015(3:4, 81:end);
naca0015_a0 = [flip(naca0015_a0_up, 2), naca0015_a0_low];
naca0015_a5 = [flip(naca0015_a5_up, 2), naca0015_a5_low];
x = linspace(0, 1, 100); % vector of an airfoil's x coordinates to be used with ThinAirfoil.m and camberline.m
%% Plotting pressure distribution about the airfoil at 0 and 5 dgr AoA
figure
h1 = axes;
hold on
plot(naca0015_a0(1, :), naca0015_a0(2, :))
plot([0, 1], [0, 0], 'k--')
set(h1, 'Ydir', 'reverse')
xlabel('x/c')
ylabel('Cp')
grid on
grid minor
hold off

dcp_naca0015_a5 = naca0015_a5_low(2, :) - naca0015_a5_up(2, :);
figure
h2 = axes;
hold on
plot(naca0015_a5(1, :), naca0015_a5(2, :))
plot(naca0015_a5_up(1, :), -dcp_naca0015_a5)
plot([0, 1], [0, 0], 'k--')
xlabel('x/c')
ylabel('Cp')
set(h2, 'Ydir', 'reverse')
text(0.05, -2.2, '\leftarrow-dCp', 'Color','r')
grid on
grid minor
hold off

%% Comparison of pressure distributions 
% obtained with thin airfoil theory and one obtained from a potential flow 
% analysis based on a surface singularity method
figure
h3 = axes;
hold on
plot(naca0015_a5_up(1, :), -dcp_naca0015_a5)
Cp_0015_thin_airfoil_mid = ThinAirfoilMidpoints(5, [x; zeros(1, length(x))]);
plot(x, -Cp_0015_thin_airfoil_mid)
legend('Surface singularity method', 'Thin airfoil theory')
xlabel('x/c')
ylabel('-dCp')
set(h3, 'Ydir', 'reverse')
ylim([-3.5, 0])
grid on
grid minor
hold off

%% Difference in upper and lower surface pressure distribution of a cambered airfoil 
z_2412 = camberline(x, 2, 4);
Cp_2412_thin_airfoil_mid = ThinAirfoilMidpoints(5, [x; z_2412]);
figure
h4 = axes;
hold on
plot(x, -Cp_2412_thin_airfoil_mid, 'b:')
plot([0, 1], [0, 0], 'k--')
xlabel('x/c')
ylabel('-dCp')
set(h4, 'Ydir', 'reverse')
grid on
grid minor
hold off

%% Lift coefficient (vs. AoA) comparison between a symmetrical and a cambered airfoils
Cl_a_0015_thin_airfoil = zeros(1, 20);
Cl_a_2412_thin_airfoil = zeros(1, 20);
for i = 1:21
    [~,~,Cl_a_0015_thin_airfoil(i)] =  ThinAirfoilMidpoints(i-11, [x; zeros(1, length(x))]);
    [~,~,Cl_a_2412_thin_airfoil(i)] = ThinAirfoilMidpoints(i-11, [x; z_2412]);
end
figure
hold on
plot(linspace(-10, 10, 21), Cl_a_0015_thin_airfoil)
plot(linspace(-10, 10, 21), Cl_a_2412_thin_airfoil)
plot([-20 20], [0 0], 'k', [0 0], [-1.5 1.5], 'k')
ylim([-1.5, 1.5])
xlim([-20, 20])
xlabel('Angle of attack [dgr]')
ylabel('Lift coefficient')
legend('Camber 0%', 'Camber 2%', 'Location', 'northwest')
grid on
grid minor
hold off


