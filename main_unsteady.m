%RWA - assignment 3 - unsteady case
clc
clear all
close all

chord = 1;
coord(:,1) = linspace(0,chord,100);  % flat plate coordinates (x,z)
coord(:,2) = zeros(1,100);
Uinf = 10;  % freestream velocity
N = 20;
rho = 1.225;

t = [0:0.1:5];
k = 0.02; %[0.02,0.05,0.1] reduced frequency
Omega = k*2*Uinf/chord;
pitch = 0.2*sin(t);% pitching motion

figure()
for i=1:length(t) % time loop
    aoa = pitch(i); % angle of attack changes (in radians)
    aoa = aoa/pi*180;
    [~,xi,zi,~,Cl(i)]=LumpedVortex(coord,chord, aoa,N, Uinf, rho,0);
    
    plot(xi,zi, 'k','LineWidth', 1.5)
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