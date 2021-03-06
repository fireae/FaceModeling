%%xx_sift can only be used under Linux system
function [desc] = local_descriptors( img, xy, dsize, dbins, options )

featType = options.descType;

stage = options.current_cascade;

dsize = options.descScale(stage) * size(img,1);

if strcmp(featType,'raw')
    
    if size(img,3) == 3
        im = im2double(rgb2gray(uint8(img)));
    else
        im = im2double(uint8(img));
    end
    npts = size(xy,2);
    for ipts = 1 : npts
        desc(ipts,:) = raw(im,xy(:,ipts), options.descRawWin*options.descRawWin);
    end

elseif strcmp(featType,'xx_sift')
    
%     i = randi([1 68],1,1);
%     rect = [xy(18,:) - [dsize/2 dsize/2] dsize dsize];
%     
%     if 1
%         figure(2); imshow(img); hold on;
%         rectangle('Position',  rect, 'EdgeColor', 'g');
%         hold off;
%         pause;
%     end
    
    
    if size(img,3) == 3
        im = im2double(rgb2gray(uint8(img)));
    else
        im = im2double(uint8(img));
    end
    
    %xy = xy - repmat(dsize/2,size(xy,1),2);
    
    desc = dSift(im,xy,dbins);
    
    %extract sift features using vl_feat
elseif strcmp(featType,'vlsift')
        if size(img,3) == 3
            im = single(rgb2gray(uint8(img)));
        else   
            im = single(img);
        end
        scale = 2;%dsize/4;
        for ipts = 1:size(xy,1)
            frame = [xy(ipts,:)';scale;0];
            [key,sift] = vl_sift(im,'frames',frame);
            desc(ipts,:) = sift';
        end
%         frame = [xy';scale*ones(1,size(xy,1));zeros(1,size(xy,1))];
%         [key,sift] = vl_sift(im,'frames',frame);
%         desc = sift';
            
    
elseif strcmp(featType,'hog')
    
    if size(img,3) == 3
        im = im2double(rgb2gray(uint8(img)));
    else
        im = im2double(uint8(img));
    end
    
    npts = size(xy,2);
    desc = zeros(npts, 2 * 2 * 31);
    for ipts = 1 : npts
        %disp(ipts);
        x = xy(1,ipts); y = xy(2,ipts);
        %if x < size(img,1) && x > 0 && y < size(img,2) && y > 0
        desc(ipts,:) = hog(im,xy(:,ipts),dsize)';
        %end
    end
    
end

end
