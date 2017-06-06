function distance = norm_Edist(gt,aligned,options)
upper_bound = min(gt(1,:));
left_bound = min(gt(2,:));
height = max(gt(1,:)) -  upper_bound;
width = max(gt(2,:)) - left_bound;

normalized_gt = [(gt(1,:) - upper_bound)/height;(gt(2,:) - left_bound)/width];
normalized_aligned = [(aligned(1,:) - upper_bound)/height;(aligned(2,:) - left_bound)/width];

tmp = sqrt(sum((normalized_gt - normalized_aligned).^2));
%tmp = sqrt(sum((gt - aligned).^2));
distance = mean(tmp);
disp(['proj dist = ' num2str(distance)]);

end