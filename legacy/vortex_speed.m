function speeds = vortex_speed(xj,zj,gamma,vertices, Uinf)
speeds = zeros(length(vertices.coord(:,1)),2);
for i=1:length(vertices.coord(:,1))
    for j=1:length(xj)  % velocities induced by airfoil
        u_ind = v_2D_from_vortex(gamma(j), [xj(j), zj(j)], vertices.coord(i,:));
        speeds(i,1) = speeds(i,1)+u_ind(1);
        speeds(i,2) = speeds(i,2)+u_ind(2);
    end
    
    for j=1:length(vertices.coord(:,1)) % velocities induced by other shed vertices
        if i~=j
            u_ind = v_2D_from_vortex(vertices.gamma(j), vertices.coord(j,:), vertices.coord(i,:));
            speeds(i,1) = speeds(i,1)+u_ind(1);
            speeds(i,2) = speeds(i,2)+u_ind(2);
        end
    end
end
speeds(:,1) = speeds(:,1)+Uinf; % should I add this???
end