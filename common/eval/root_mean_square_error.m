%%This function calculate the Vertex Distance Cost(VDC) between the ground
%%truth 3D shape and reconstructed 3D shape. It assumed that these two
%%shape has been aligned.
%%The input are [n,1] vectors with[x1;y1;z1;x2;y2;z2;...;xn;yn;zn]; 

function rmse = root_mean_square_error(gt, est, validpoint)
if(numel(gt)~=numel(est)) 
    disp('VDC:input vectors has different numbers of points!');
    pause;
end

% gt = reshape(gt, [ 3 numel(gt)/3 ])'; 
% est = reshape(est,[ 3 numel(gt)/3 ])';

%only calculate valid points 
gt = gt(:,logical(validpoint));%[3 N]
est = est(:,logical(validpoint));%[3 N]

rmse = sqrt(sum(sum((gt - est).^2))./size(gt,2));
% disp(['vertex dist = ' num2str(error)])

end