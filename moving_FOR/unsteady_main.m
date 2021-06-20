clc; clear all; close all;

Uinf = 10;
Winf = 0;
Qinf = sqrt(Uinf^2+Winf^2);
rho = 1.255;
chord = 1;
Npan = 10;  % number of panels
k = 0.05; %[0.02,0.05,0.1] reduced frequency
omega = k*2*Uinf/chord;
dt = 0.1;
Nts = 70;  % number of time steps
it = 0; % time step counter t = it*dt

% discretization
panels = plate_discretization(chord, Npan);

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

%% first time step - steady calculation to get wake vorticity
% for movement to the left at speed Uinf and pitch osciallation at frequency omega
t = 0;
kinematics.X0 = -Uinf*t;   % position of LE
kinematics.X0dot = -Uinf;  % X speed
kinematics.Z0 = 0;     % position of LE
kinematics.Z0dot = 0;  % Z speed
kinematics.theta = 0.25*sin(omega*t);
kinematics.thetadot = 0.25*omega*cos(omega*t);

vertices.gamma = [];
vertices.x = [];
vertices.z = [];   
[Cl_0,~,~,gamma_old] = LumpedVortex(1, Npan, kinematics, coeff, panels, chord, Qinf, rho, vertices);

figure()
% hold on
TE_new = [cos(kinematics.theta), sin(kinematics.theta); -sin(kinematics.theta), cos(kinematics.theta)]*[panels.x(end);panels.z(end)]+[kinematics.X0;kinematics.Z0] ;    %location of TE at t=dt
for it=1:Nts
% old TE position
TE_old = TE_new;
    
t = it*dt;
kinematics.X0 = -Uinf*t;   % position of LE
kinematics.theta = 0.25*sin(omega*t);
kinematics.thetadot = 0.25*omega*cos(omega*t);

TE_new = [cos(kinematics.theta), sin(kinematics.theta); -sin(kinematics.theta), cos(kinematics.theta)]*[panels.x(end);panels.z(end)]+[kinematics.X0;kinematics.Z0];    %location of TE at t=dt
temp = [cos(kinematics.theta), sin(kinematics.theta); -sin(kinematics.theta), cos(kinematics.theta)]*[panels.x;panels.z]+[kinematics.X0;kinematics.Z0];    %location of TE at t=dt
% clf
% hold on

[Cl(it),~,vertices, gamma] = LumpedVortex(0, Npan, kinematics, coeff, panels, chord, Qinf, rho, vertices, TE_old, TE_new, gamma_old, dt);
gamma_old = gamma;
% convect vertices - just with freestream, does not calculate induced
% velocity
plot(temp(1,:), temp(2,:),'k', 'LineWidth',1.5)
hold on
scatter(vertices.x, vertices.z,10,'filled')
hold off
axis equal
ylim([-1 1])
xlim([min(temp(1,:))-chord max(vertices.x)])
drawnow;
% hold off
pause(0.2)
vertices.x = vertices.x+Uinf*dt;
vertices.z = vertices.z+Winf*dt;
end
% hold off
figure()
plot(1:Nts,rad2deg(0.25*sin([1:Nts]*omega*dt)))
figure()
plot(1:Nts, Cl)
figure()
plot(rad2deg(sin([1:Nts]*omega*dt)),Cl)
