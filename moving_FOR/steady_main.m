clc; clear all; close all;
Uinf = 10;
Winf = 0;
Qinf = sqrt(Uinf^2+Winf^2);
rho = 1.255;
t = 0;
omega = 0;

% discretization
Npan = 10;
chord = 1;
panels = plate_discretization(chord, Npan);

% for movement to the left at speed Uinf and pitch osciallation at
% frequency omega
kinematics.X0 = -Uinf*t;   % position of LE
kinematics.X0dot = -Uinf;  % X speed
kinematics.Z0 = 0;     % position of LE
kinematics.Z0dot = 0;  % Z speed
kinematics.thetadot = 0;

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

Cl = [];
for theta=0%deg2rad([-15:1:15])
    kinematics.theta = theta;
%     temp = [cos(kinematics.theta), sin(kinematics.theta); -sin(kinematics.theta), cos(kinematics.theta)]*[panels.x;panels.z]+[kinematics.X0;kinematics.Z0] ;    %location of TE at t=dt
%     X = temp(1,:);
%     Z = temp(2,:);
    
    [temp] = LumpedVortex(1, Npan, kinematics, coeff, panels, chord, Qinf, rho);
    Cl = [Cl,temp];
end

plot([-15:1:15], Cl)
grid on
grid minor
xlabel("\alpha [^\circ]")
ylabel("C_L")
