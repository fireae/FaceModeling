function projection = cal_weak_perspective( point3d,s,R,T)
%CAL_WEAK_PERSPECTIVE Summary of this function goes here
%   Detailed explanation goes here

A = [1 0 0;0 1 0];
projection = (A*(s.*R*point3d +repmat(T,1,size(point3d,2))));


end

