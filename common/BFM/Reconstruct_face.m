function Rv = Reconstruct_face(model, para, trans)

Re_shape  = coef2object(para, model.shapeMU, model.shapePC, model.shapeEV );
Re_shape = reshape(Re_shape,3,[]);
R = Ang2Rot(trans(1:3));  
 
 Trans(1:3,1:3) = R;
 Trans(1:3,4)    = trans(4:6);
 Trans(4,:) = [0,0,0,1];
 
 Re_shape = Trans * [Re_shape;ones(1,size(Re_shape,2))];
 Rv = Re_shape(1:3,:);
 
%  Rn = zeros(size(Rv));
%  faces = model.tl;
%  tmp = Rn(:,faces(:,1));
%  Rn(:,faces(:,1)) = tmp + cross(Rv(:,faces(:,3)) - Rv(:,faces(:,1)),Rv(:,faces(:,2)) - Rv(:,faces(:,1)));
%  
%  test = cross(Rv(:,faces(:,3)) - Rv(:,faces(:,1)),Rv(:,faces(:,2)) - Rv(:,faces(:,1)));
%  faces1 = faces(:,1);
%  index = 3;
%  ans =  find(faces1 == index);
%  result = sum(test(:,ans),2);
% 
%  
%  Rn(:,faces(:,2)) = Rn(:,faces(:,2)) + cross(Rv(:,faces(:,1)) - Rv(:,faces(:,2)),Rv(:,faces(:,3)) - Rv(:,faces(:,2)));
%  Rn(:,faces(:,3)) = Rn(:,faces(:,3)) + cross(Rv(:,faces(:,2)) - Rv(:,faces(:,3)),Rv(:,faces(:,1)) - Rv(:,faces(:,3)));
% 

end