function panels = plate_discretization(chord, N, flap, beta, flap_c)
% creates flat plate
panels.x = 0:1/N:chord;
panels.z = zeros(1,N+1);
if flap & nargin>3
    [~,id] = find(panels.x >=(1-flap_c)*chord); % find the hinge
    panels.z(id) = -sin(beta)*(panels.x(id)-(1-flap_c)*chord);% for all panels behind the hinge
    panels.x(id) = (1-flap_c)*chord+cos(beta)*(panels.x(id)-(1-flap_c)*chord);% for all panels behind the hinge
end

panels.xj = zeros(1, N);
panels.xi = zeros(1, N);
panels.zj = zeros(1, N);
panels.zi = zeros(1, N);
panels.a = zeros(1, N);
panels.n = zeros(N, 2);
for i =1:N      % loop over each panel
    % Vortex points
    panels.xj(i)= panels.x(i)+1/4*(panels.x(i+1)-panels.x(i));
    panels.zj(i)= panels.z(i)+1/4*(panels.z(i+1)-panels.z(i));
    
    % Collocation points
    panels.xi(i)= panels.x(i)+3/4*(panels.x(i+1)-panels.x(i));
    panels.zi(i)= panels.z(i)+3/4*(panels.z(i+1)-panels.z(i));
        
    % Normal vectors and panel slopes in collocation points
    panels.a(i)=atan2((panels.z(i+1)-panels.z(i)),(panels.x(i+1)-panels.x(i)));    % angle of panel used for calculating normal vector
    panels.n(i,:) = [cos(panels.a(i)), sin(panels.a(i))];
    panels.n(i,:)= ([cos(pi/2), -sin(pi/2); sin(pi/2), cos(pi/2)]*panels.n(i,:)')';   %normal vector of i-th panel
    panels.a(i)=asin(panels.n(i,1));  % angle of panel consistent with the following equations
end

panels.tau = [cos(panels.a'), -sin(panels.a')]; % tangential components

panels.eta = -(panels.zi+panels.zj);    % camber - at middle of panels, so average of collocation and vortex points
panels.deta = zeros(1,N);   % time change of camber
end