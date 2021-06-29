%% checking convergence - reaction to step change
clc; clear all; close all;

Uinf = 10;
Winf = 0;
Qinf = sqrt(Uinf^2+Winf^2);
rho = 1.255;
chord = 1;
Npan = 10;  % number of panels
k = 0.1; %[0.02,0.05,0.1] reduced frequency
omega = k*2*Uinf/chord;
theta0 = deg2rad(10); % here - step change in angle of attack
dt = 0.1;
Nts = 100;  % number of time steps
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
t = dt;
vertices.gamma = [];
vertices.x = [];
vertices.z = [];
gamma_old = zeros(Npan,1);
p_atm = 0; %101300;
Ngrid = 100;
TE_new = [1,0];
% TE_new = [cos(kinematics.theta), sin(kinematics.theta); -sin(kinematics.theta), cos(kinematics.theta)]*[panels.x(end);panels.z(end)]+[kinematics.X0;kinematics.Z0] ;    %location of TE at t=dt
% kinematics.thetadot = 0;

for it=1:Nts
    kinematics.X0 = -Uinf*t;   % position of LE
    kinematics.X0dot = -Uinf;  % X speed
    kinematics.Z0 = 0;     % position of LE
    kinematics.Z0dot = 0;  % Z speed
    kinematics.theta = theta0;
    kinematics.thetadot = 0;%theta0/dt;
    % old TE position
    TE_old = TE_new;
    
    t = (it-1)*dt;
    kinematics.X0 = -Uinf*t;   % position of LE
%     kinematics.theta = theta0*sin(omega*t);
%     kinematics.thetadot = theta0*omega*cos(omega*t);
    
    TE_new = [cos(kinematics.theta), sin(kinematics.theta); -sin(kinematics.theta), cos(kinematics.theta)]*[panels.x(end);panels.z(end)]+[kinematics.X0;kinematics.Z0];    %location of TE at t=dt
    temp = [cos(kinematics.theta), sin(kinematics.theta); -sin(kinematics.theta), cos(kinematics.theta)]*[panels.x;panels.z]+[kinematics.X0;kinematics.Z0];    %location of TE at t=dt
    
    
    [Cl(it),~,vertices, gamma, Cl_c(it), Cl_nc(it)] = LumpedVortex(0, Npan, kinematics, coeff, panels, chord, Qinf, rho, vertices, TE_old, TE_new, gamma_old, dt);
    gamma_old = gamma(1:end-1);
    plate_vert.gamma = [gamma(1:end-1);vertices.gamma'];   % vertices on the panels
    plate_vert.coord = [[cos(kinematics.theta), sin(kinematics.theta); -sin(kinematics.theta), cos(kinematics.theta)]*[panels.xj;panels.zj]+...
        [kinematics.X0;kinematics.Z0], [vertices.x;vertices.z]] ;
    
    % convect vertices - just with freestream, does not calculate induced velocity
    vertices.x = vertices.x+Uinf*dt;
    vertices.z = vertices.z+Winf*dt;
end

figure()
hold on
plot([1:Nts]*dt, Cl_c)
plot([1:Nts]*dt, Cl)
hold off
grid on
grid minor
xlabel("t")
ylabel("C_L")

