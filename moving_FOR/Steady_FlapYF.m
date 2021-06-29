%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%             START          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

beta=30;
N_mesh=100;                      % Number of points on the mesh
N_foil=30;                       % Number of points on the airfoil
n_foil=N_foil-1;                 % Number of sections
n_mesh=N_mesh-1;
alpha=10;
CL_current=[];
CL_theory=[];
CM=[];
styles = {'-', '--', '-.',':', '-','--'};
st=1;
figure()
%plot([0:0.5:10],2*pi*sind([0:0.5:10]),'-r+')
for beta=[0 5 10 20 30]
  CL_current=[];
  for alpha=[-10:1:10]
    alpha    
    rho=1.225;
    U_inf=10;
    c=1;                         % Chord length
    x=linspace(0,3*c,N_mesh);
    z=linspace(0,3*c,N_mesh);
    [X,Z]=meshgrid(x,z);         % Generating X and Y matrices

    %% Gathering airfoil coordinates
    x_c4=1*c;
    z_c4=1.5*c;

    mat_rot=[cosd(alpha) -sind(alpha) ; sind(alpha) cosd(alpha)];
    [x_vp,z_vp,x_cp,z_cp,LE_vec,TE_vec]=airfoil_geometry(c,x_c4,z_c4,mat_rot,N_foil,n_foil,beta);
    d0=f_d0(c,x_c4,z_c4,mat_rot,N_foil,n_foil);


    %% Set the normals of the 
    dx_local=abs(x_vp-x_cp);
    dz_local=z_vp-z_cp;
    alpha_local=atand(dz_local./dx_local);
    normals=[sind(alpha_local)' cosd(alpha_local)'];

    %% Set up influence matrix
    inf_mat=zeros(n_foil,n_foil);
    gamma_set=1;
    for i=[1:n_foil]
        for j=[1:n_foil]
            [u_ind,w_ind] = VOR2D(x_cp,x_vp,z_cp,z_vp,i,j,gamma_set);
            inf_mat(i,j)=u_ind*normals(i,1)+w_ind*normals(i,2);
        end 
    end
    RHS=-U_inf*normals(:,1);
    Gamma=inf_mat\RHS;
    dL=rho*U_inf*Gamma;
    CL_current=[CL_current sum(dL)/(0.5*rho*U_inf^2*c)];
    CL_theory=[CL_theory 2*pi*sind(alpha)];
    inner=dL.*d0';
    CM=[CM -sum(inner)/(0.5*rho*U_inf^2*c^2)];

    %% Determine velocity and pressure field
    u_field=zeros(N_mesh^2,1);
    w_field=zeros(N_mesh^2,1);
    for i=[1:N_mesh^2]
        for j=[1:n_foil]
            [u_ind_point,w_ind_point] = VOR2D(X(:)',x_vp,Z(:)',z_vp,i,j,Gamma(j));
            u_field(i)=u_field(i)+u_ind_point;
            w_field(i)=w_field(i)+w_ind_point;
        end 
    end 
    u_field=reshape(u_field,N_mesh,N_mesh)+U_inf;
    w_field=reshape(w_field,N_mesh,N_mesh);

    V_field=sqrt(u_field.^2+w_field.^2);
    p_field=0.5*V_field.^2*rho;
    p_field=0.5*rho*U_inf^2-p_field;
  end 
  hold on
  plot([-10:1:10],CL_current, 'linestyle', styles{st})
  st=st+1;
end
%%
q=legend('\beta=0^\circ','\beta=5^\circ','\beta=10^\circ','\beta=20^\circ','\beta=30^\circ');

% set(q,'Interpreter','LaTex','location','northwest')
grid
xlabel('\alpha [^\circ]') %,'fontsize',22,'Interpreter','LaTex')
ylabel('C_L [-]')%,'fontsize',22,'Interpreter','LaTex')

figure()
hold on
h = pcolor(X,Z,V_field);
set(h, 'EdgeColor', 'none');
colormap pink(32)
g=colorbar;
caxis([3 20])
shading interp
title("Velocity field")
xlabel("x")
ylabel("z")
% xlabel(g,'Velocity [m/s]','fontsize',22,'Interpreter','LaTex')
xlim([0,3])
ylim([0,3])
plot([LE_vec(1) x_cp TE_vec(1)],[LE_vec(2) z_cp TE_vec(2)],'k-','Linewidth',4)
t=16;
quiver(X(1:t:end),Z(1:t:end),u_field(1:t:end),w_field(1:t:end),1,'color',[0 0 0])

figure()
hold on
h = pcolor(X,Z,p_field);
set(h, 'EdgeColor', 'none');
colormap pink(32)
g=colorbar;
caxis([-150 100])
shading interp
title("Pressure field")
xlabel("x")
ylabel("z")
% xlabel(g,'Pressure [Pa]','fontsize',22,'Interpreter','LaTex')
xlim([0,3])
ylim([0,3])
plot([LE_vec(1) x_cp TE_vec(1)],[LE_vec(2) z_cp TE_vec(2)],'k-','Linewidth',4)




















