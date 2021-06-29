function [x_vp,z_vp,x_cp,z_cp,LE_vec,TE_vec]=airfoil_geometry(c,x_c4,z_c4,mat_rot,N_foil,n_foil,beta)
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

%Set flap
flap_index=find(x_vp>=0.5*c);
distance_rot_point_vp=x_vp(flap_index)-0.5*c;
distance_rot_point_cp=x_cp(flap_index)-0.5*c;

x_vp(flap_index)=0.5*c+cosd(beta)*distance_rot_point_vp;
x_cp(flap_index)=0.5*c+cosd(beta)*distance_rot_point_cp;

z_vp(flap_index)=sind(beta)*distance_rot_point_vp;
z_cp(flap_index)=sind(beta)*distance_rot_point_cp;

TE_dist=x_TE-0.5*c;
x_TE=0.5*c+TE_dist*cosd(beta);
z_TE=sind(beta)*TE_dist;

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

LE_vec=mat_rot*[x_LE ; z_LE];
LE_vec(1)=LE_vec(1)+x_c4;
LE_vec(2)=-LE_vec(2)+z_c4;

x_vp=x_vp+x_c4;
z_vp=z_vp+z_c4;
x_cp=x_cp+x_c4;
z_cp=z_cp+z_c4;
end 
