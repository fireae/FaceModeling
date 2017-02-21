%%This function calculate the Vertex Distance Cost(VDC) between the ground
%%truth 3D shape and reconstructed 3D shape. It assumed that these two
%%shape has been aligned.
%%The input are [n,1] vectors with[x1;y1;z1;x2;y2;z2;...;xn;yn;zn]; 

function error = vdc(ground_truth, aligned_shape)
if(numel(ground_truth)~=numel(aligned_shape)) 
    disp('VDC:input vectors has different numbers of points!');
    pause;
end

ground_truth = reshape(ground_truth, [ 3 numel(ground_truth)/3 ])'; 
aligned_shape = reshape(aligned_shape,[ 3 numel(ground_truth)/3 ])';

% min_p = min(ground_truth,[],1);
% max_p = max(ground_truth,[],1);

%%normalize ground truth ans aligned shape
% normalized_gt(:,1) = (ground_truth(:,1) - min_p(1))/(max_p(1) - min_p(1));
% normalized_gt(:,2) = (ground_truth(:,2) - min_p(2))/(max_p(2) - min_p(2));
% normalized_gt(:,3) = (ground_truth(:,3) - min_p(3))/(max_p(3) - min_p(3));
% 
% normalized_as(:,1) = (aligned_shape(:,1) - min_p(1))/(max_p(1) - min_p(1));
% normalized_as(:,2) = (aligned_shape(:,2) - min_p(2))/(max_p(2) - min_p(2));
% normalized_as(:,3) = (aligned_shape(:,3) - min_p(3))/(max_p(3) - min_p(3));
% 

%error = mean(sqrt(sum((normalized_gt - normalized_as).^2,2)));
error = mean(sqrt(sum((ground_truth - aligned_shape).^2,2)));

end