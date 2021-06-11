%RWA - assignment 3 - unsteady case
%%CHANGE!!! PITCHING MOVEMENT IS NOT JUST THE ANGLE! 
clc
clear all
close all

chord = 1;
coord(:,1) = linspace(0,chord,100);  % flat plate coordinates (x,z)
coord(:,2) = zeros(1,100);
Uinf = 10;  % freestream velocity
N = 20;
rho = 1.225;

t = [0:1:10];
k = 0.02; %[0.02,0.05,0.1] reduced frequency
Omega = k*2*Uinf/chord;
theta0 = 1;
pitch = theta0*sin(Omega*t);% pitching motion


% first calculation for t=0
aoa = pitch(1); % angle of attack changes (in radians)
aoa = aoa/pi*180;
[~,xi,zi,~,Cl(1), circ_old, vortex_coord]=LumpedVortex(coord,chord, aoa,N, Uinf, rho,0,1,0, [0,0], pitch(1));

figure()
hold on
plot(xi,zi, 'k','LineWidth', 1.5)
scatter(vortex_coord(1), vortex_coord(2))
hold off
xlim([-chord 2*chord])
ylim([-chord chord])
xlabel("x")
ylabel("z")
getframe(gcf);


for i=2:length(t) % time loop
    aoa = pitch(i); % angle of attack changes (in radians)
    aoa = aoa/pi*180;
    [~,xi,zi,~,Cl(i),circ_old, vortex_coord]=LumpedVortex(coord,chord, aoa,N, Uinf, rho,0,1,circ_old,vortex_coord, pitch(i)-pitch(i-1));
    
    hold on 
    plot(xi,zi, 'k','LineWidth', 1.5)
    scatter(vortex_coord(1), vortex_coord(2))
    hold off
    xlim([-chord 2*chord])
    ylim([-chord chord])
    xlabel("x")
    ylabel("z")
    getframe(gcf);
end

figure()
plot(t,Cl)
grid on
grid minor
xlabel("time [s]")
ylabel ("Lift coefficient C_L")

figure()
plot(pitch/pi*180,Cl)
grid on
grid minor
xlabel("Angle of attack [^\circ]")
ylabel ("Lift coefficient C_L")