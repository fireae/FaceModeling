function point2D = cal_key_projection(Points3N, cx, cy, fx, fy, Rc, Tc)
 %Calcluate the projection positions of 3D face model keypoints according
 %to the parameters, which include 200 shape parameters, 3 rotation paras,
 %3 translation paras, focal length
  
  TransVi = Rc * Points3N +  repmat(Tc,[1,size(Points3N,2)]);
  point2D(:,1) = round(cx - fx*TransVi(1,:)./TransVi(3,:));
  point2D(:,2) =round(cy - fy*TransVi(2,:)./TransVi(3,:));
 %keyIndx = reshape([3*keyIndx+ones(size(keyIndx));3*keyIndx+2*ones(size(keyIndx));3*keyIndx+3*ones(size(keyIndx))],[],1);
%  Mi = model.shapeMU(keyIndx,:);
%  Pi = model.shapePC(keyIndx,:);
%  Vi = Mi + Pi*para(1:199);
%  Vi = reshape(Vi,3,[]);
% %  tmp(1:2,:) = Vi(2:3,:);
% %  tmp(3,:) = Vi(1,:);
% %  Vi = tmp;
% %  clear tmp;
%  R = Ang2Rot(para(201:203));  
%  
%  Trans(1:3,1:3) = R;
%  Trans(1:3,4)    = para(204:206);
%  Trans(4,:) = [0,0,0,1];
%  
%  Vi = Trans * [Vi;ones(1,size(Vi,2))];
%  TransVi = Vi(1:3,:)./[Vi(4,:);Vi(4,:);Vi(4,:)];

  
  
end  