function desc = hog( im, pos , lmsize )

%fsize  = sqrt(norm_size);
%lmsize  = fsize;
%gsize = options.canvasSize(1) * options.descScale(1);

rect =  [pos(1) - (lmsize-1)/2, ...
         pos(2) - (lmsize-1)/2, ...
         lmsize - 1, lmsize - 1];
     


cropim = imcrop(im,rect);
if 0
   figure; imshow(cropim);hold on;
   pause;
end
%disp([size(cropim) lmsize]);
scalesize = 64;
%scale image patch to [64 64]
if size(cropim,1) ~= scalesize || size(cropim,2) ~=scalesize
     cropim = imresize(cropim,[scalesize scalesize]);
end

cellSize = 32 ;
%tmp = vl_hog(single(cropim), cellSize, 'verbose');
tmp = vl_hog(single(cropim), cellSize);

if 0
   figure; imshow(tmp);
   pause;
end

%desc = feat_normalize(tmp(:));
desc = tmp(:);
     
end