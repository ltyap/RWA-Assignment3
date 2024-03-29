clc; clear all; close all;

Uinf = 1;
Winf = 0;
Qinf = sqrt(Uinf^2+Winf^2);
rho = 1.255;
chord = 1;
Npan = 10;  % number of panels
k = [0.02,0.05,0.1]; %reduced frequency
omega = k*2*Uinf/chord;
theta0 = 0.175;
dt = 0.1;
Nts = 200;  % number of time steps
t = [0:Nts]*dt;
% theta = theta0*sin(omega*t);
styles = {'-', '--', '-.',':', '-','--'};
colors = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.4940, 0.1840, 0.5560];[0.4660, 0.6740, 0.1880]]; % define color scheme, more colors need to be added for larger analysis

figure()
hold on
for i=1:length(k)
    plot(t,theta0*sin(omega(i)*t),'linestyle',styles{i}, 'Color', colors(i,:), 'DisplayName',sprintf('k =%.2f', k(i)))
end
hold off
grid on
grid minor
legend show
xlabel("t")
ylabel("\theta")

%% different k
clear t
for i=1:length(k)
    temp= load("unsteady_cl_theta_k"+k(i)+".mat");
    t(i,:) = temp.Cl(1,:);
% %     temp= load("unsteady_theta_k"+k(i)+".mat");
    theta(i,:) = theta0*sin(k(i)*2*Uinf/chord*t(i,:));
%     temp= load("unsteady_cl_theta_k"+k(i)+".mat");
    cl(i,:) = temp.Cl(2,:);
end

figure()
hold on
for i=1:length(k)
    plot(theta(i,:),cl(i,:),'linestyle',styles{i}, 'Color', colors(i,:), 'DisplayName',sprintf('k =%.2f', k(i)))
%     h = quiver(theta(i,1),cl(i,1), theta(i,2)-theta(i,1), cl(i,2)-cl(i,1),1.5,'filled', 'Color', colors(i,:))
%     h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
hold off
grid on
grid minor
legend show
xlabel("\theta")
ylabel("C_L")

s=2*t*Uinf/chord; % in distance semichords
figure()
hold on
for i=1:length(k)
    plot(s(i,:),cl(i,:),'linestyle',styles{i}, 'Color', colors(i,:), 'DisplayName',sprintf('k =%.2f', k(i)))
end
hold off
grid on
grid minor
legend show
xlim ([0,250])
xlabel("reduced time")
ylabel("C_L")
