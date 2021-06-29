function [d0]=f_d0(c,x_c4,z_c4,mat_rot,N_foil,n_foil)
x_LE=-c/4;
x_TE=3*c/4;
z_LE=0;
z_TE=0;

% Airfoil discretisation
x_dp=linspace(x_LE,x_TE,N_foil);           % dp=discretisation point
z_dp=linspace(z_LE,z_TE,N_foil);

% Set control point and vortex points
element_length=c/(N_foil-1);
x_vp=x_dp(1:end-1)+element_length/4;
x_cp=x_dp(1:end-1)+3*element_length/4;
z_vp=z_dp(1:end-1);
z_cp=z_dp(1:end-1);


d0=sqrt(x_cp.^2+z_cp.^2);
index=find(x_cp<0);
d0(index)=-d0(index);

for i=[1:n_foil]
    vp_vec=mat_rot*[x_vp(i);z_vp(i)];
    cp_vec=mat_rot*[x_cp(i);z_cp(i)];
    
    x_vp(i)=vp_vec(1);
    z_vp(i)=-vp_vec(2);
    x_cp(i)=cp_vec(1);
    z_cp(i)=-cp_vec(2);
    
    
end 
TE_vec=mat_rot*[x_TE ; z_TE];
TE_vec(1)=TE_vec(1)+x_c4;
TE_vec(2)=-TE_vec(2)+z_c4;

x_vp=x_vp+x_c4;
z_vp=z_vp+z_c4;
x_cp=x_cp+x_c4;
z_cp=z_cp+z_c4;
end 
