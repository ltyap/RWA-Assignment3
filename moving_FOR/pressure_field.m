function p_grid = pressure_field(Uinf, p_atm, rho, x_grid, z_grid, u_grid, v_grid)
p_grid = zeros(size(x_grid));
% from bernoulli equation - for a steady case!!!
dp = rho*((u_grid.^2+v_grid.^2)/2 - Uinf^2/2);
p_grid = p_atm-dp;
end