% RWA - assignment 3
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
    [~,~, Cl(i)]=LumpedVortex(coord,chord, aoa(i),N, Uinf, rho,0);
end

figure()
plot(aoa,Cl)
grid on
grid minor
xlabel("Angle of attack [^\circ]")
ylabel ("Lift coefficient C_L")
