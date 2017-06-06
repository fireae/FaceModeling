%%This function calculate the Vertex Distance Cost(VDC) between the ground
%%truth 3D shape and reconstructed 3D shape. It assumed that these two
%%shape has been aligned.
%%The input are [n,1] vectors with[x1;y1;z1;x2;y2;z2;...;xn;yn;zn]; 

function error = euclidian_dist(ground_truth, aligned_shape)
if(numel(ground_truth)~=numel(aligned_shape)) 
    disp('VDC:input vectors has different numbers of points!');
    pause;
end

error = sqrt(sum((ground_truth - aligned_shape).^2));
% disp(['vertex dist = ' num2str(error)])

end