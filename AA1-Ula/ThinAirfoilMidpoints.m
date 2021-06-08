function [dCp, gamma, Cl] = ThinAirfoilMidpoints(alpha, coord, V, rho)
    alpha = alpha*pi/180;
    x_pre = coord(1, :);
    z_pre = coord(2, :);
    x=zeros(1,length(x_pre));
    z=zeros(1,length(x_pre));
    
    for i=1:length(x_pre)-1
     	x(i) = x_pre(i)+0.5*(x_pre(i+1)-x_pre(i));
        z(i)=z_pre(i)+0.5*(z_pre(i+1)-z_pre(i));
    end
    x(length(x_pre)) = 1;
    z(length(x_pre)) = 0;
    
    theta = acos(1-2*x);
    dCp = zeros(1, length(theta));
    dzdx = zeros(1, length(theta));
    gamma = zeros(1, length(theta));
    for i = 1:length(theta)-1
        dzdx(i) = (z(i+1)-z(i))/(x(i+1)-x(i));
    end
    A0 = alpha - 1/pi * trapz(theta, dzdx);
    A = zeros(5, 1);
    for i = 1:length(A)
        A(i) = 2/pi * trapz(theta, dzdx.*cos(i*theta));
    end
    for i = 1:length(theta)
        sum = 0;
        for j = 1:length(A)
            sum = sum + (A(j)*sin(j*theta(i)));
        end
        gamma(i) = 2*V*(A0*(1+cos(theta(i)))/sin(theta(i)) + sum);
        dCp(i) = rho*V*gamma(i)/(0.5*rho*V^2);
    end
    alpha0 = -1/pi*trapz(theta, dzdx.*(cos(theta)));
    Cl = 2*pi*(alpha-alpha0);
end