%% RWA - assignment 3 - unsteady case
clc
clear all
close all
%% 
chord = 1;
coord(:,1) = linspace(0,chord,100);  % flat plate coordinates (x,z)
coord(:,2) = zeros(1,100);
Uinf = 10;  % freestream velocity
Npan = 40;
rho = 1.225;

x = [0,chord/4,chord];
z = [0,0,0];

dt = 0.1;
t = [0:dt:5];
k = 0.05; %[0.02,0.05,0.1] reduced frequency
Omega = k*2*Uinf/chord;
theta0 = 0.25;
pitch = theta0*sin(Omega*t);% pitching motion

% first time step
aoa = pitch(1);
[p,Cl(1),vertices]=LumpedVortex_different(coord, chord, aoa, Npan, Uinf, rho,0, 1, pitch(1));
[p_st,Cl_st(1),vertices_st]=LumpedVortex_different(coord, chord, aoa, Npan, Uinf, rho,0);


% convect vertices according to their speeds
vertices.coord(:,1) = vertices.coord(:,1)+Uinf*dt;%vertices.speed(:,1)*dt; % x coordinate
% vertices.coord(:,2) = vertices.coord(:,2)+vertices.speed(:,2)*dt; % z coordinate

zi = (chord/4-x)*sin(aoa);
for i=1:length(x)
    if x(i)<chord/4
        xi(i) = chord/4-(chord/4-x(i))*cos(aoa);
    else
        xi(i) = chord/4+(x(i)-chord/4)*cos(aoa);
    end
end

figure()
plot(xi,zi, 'k','LineWidth', 1.5)
xlim([-chord 2*chord])
ylim([-chord chord])
xlabel("x")
ylabel("z")
getframe(gcf);


for i=2:length(t) % time loop
    aoa = pitch(i); % angle of attack changes (in radians)
    [p,Cl(i),vertices]=LumpedVortex_different(coord,chord, aoa,Npan, Uinf, rho,0,1,pitch(i)-pitch(i-1),vertices);
    
    % convect vertices according to their speeds (NO - simplification - convects only with Uinf)
    vertices.coord(:,1) = vertices.coord(:,1)+Uinf*dt; %vertices.speed(:,1)*dt; % x coordinate
%     vertices.coord(:,2) = vertices.coord(:,2)+vertices.speed(:,2)*dt; % z coordinate
    
    zi = (chord/4-x)*sin(aoa);
    for i=1:length(x)
        if x(i)<chord/4
            xi(i) = chord/4-(chord/4-x(i))*cos(aoa);
        else
            xi(i) = chord/4+(x(i)-chord/4)*cos(aoa);
        end
    end
    
    hold on
    plot(xi,zi, 'k','LineWidth', 1.5)
    scatter(vertices.coord(:,1), vertices.coord(:,2))
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
plot(rad2deg(pitch),Cl)
grid on
grid minor
xlabel("Angle of attack [^\circ]")
ylabel ("Lift coefficient C_L")