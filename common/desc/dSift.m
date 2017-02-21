function desc = dSift(im, xy, dbins)
pos = [min(xy(:,1)), min(xy(:,2))];
w = max(xy(:,1)) - min(xy(:,1));
h = max(xy(:,2)) - min(xy(:,2));
rect = [pos(1) - 2*dbins, pos(2) - 2*dbins, h+4*dbins, w+4*dbins];
%rect = [pos(1), pos(2) , w, h];
% figure(2);
% imshow(im);
% hold on;plot(xy(:,1),xy(:,2),'.r');
% hold off;
cropim = imcrop(im,rect);
cropim = vl_imsmooth(cropim, sqrt((dbins/3)^2 - .25)) ;
[idx, denseSift] = vl_dsift(single(cropim), 'size', dbins,'fast','FloatDescriptors');
% figure(1)
% imshow(cropim);
% hold on;
%plot(xy(:,1)-(pos(1)-1)*ones(size(xy,1),1),xy(:,2)-(pos(2)-1)*ones(size(xy,1),1),'.r');
% hold on;
% plot(idx(1,:),idx(2,:),'r.');
% plot(xy(:,1)-(pos(1)-8)*ones(size(xy,1),1),xy(:,2)-(pos(2)-8)*ones(size(xy,1),1),'.g');
% hold off;
index = (xy(:,2) - (pos(2)-2) * ones(size(xy,1),1)) * (w+4) + xy(:,1) - (pos(1)-3) * ones(size(xy,1),1) +ones(size(xy,1),1);
desc = denseSift(:, index);
%close all;
end