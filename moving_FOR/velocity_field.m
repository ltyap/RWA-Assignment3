function [u_grid, w_grid] = velocity_field(Uinf, Winf, x_grid,z_grid,vertices)
u_grid = zeros(size(x_grid));
w_grid = zeros(size(x_grid));

for i=1:length(x_grid(:,1))
    for j=1:length(x_grid(1,:)) % loop over all grid points
        vel_ind=zeros(2,1);
        for k=1:length(vertices.gamma)    % loop over all vortex points
            temp = v_2D_from_vortex(vertices.gamma(k), vertices.coord(:,k)', [x_grid(i,j), z_grid(i,j)]);
            vel_ind = vel_ind+temp;
        end
        u_grid(i,j) = vel_ind(1);
        w_grid(i,j) = vel_ind(2);
    end
end

u_grid = u_grid + Uinf;
w_grid = w_grid + Winf;
end