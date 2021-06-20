function [xi,zi,xj,zj,n] = rotate_aoa_vector(xi,zi,xj,zj,c,aoa)
% rotates flat plate by a given angle of attack, wrt c/4
% aoa in radians

zi = (c/4-xi)*sin(aoa);
zj = (c/4-xj)*sin(aoa);
for i=1:length(xi)
    if xi(i)<c/4
    xi(i) = c/4-(c/4-xi(i))*cos(aoa);
    else
    xi(i) = c/4+(xi(i)-c/4)*cos(aoa);
    end
    
    if xj(i)<c/4
    xj(i) = c/4-(c/4-xj(i))*cos(aoa);
    else
    xj(i) = c/4+(xj(i)-c/4)*cos(aoa);
    end
end

% rotate normals
n(:,1) = sin(aoa)*ones(length(xi),1);
n(:,2) = cos(aoa)*ones(length(xi),1);

end