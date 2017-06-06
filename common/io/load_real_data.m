function Data = load_real_data(path, options)
%% output format
%{
DATA.
- width_orig: the width of the original image.
- height_orig: the height of the original image.
- img_gray: the crop image.
- height: the height of crop image.
- wdith: the width of crop image.
- shape_gt: ground-truth landmark.
- bbox_gt: bounding box of ground-truth.
- bbox_facedet: face detection region
%}

imlist = dir([path '*_N_*png']);
nimgs  = length(imlist);
%nimgs  = 10;

%% face detection
% Create a cascade detector object.
faceDetector = vision.CascadeObjectDetector();
%faceDetector.MergeThreshold = 7;
Data = cell(nimgs, 1);

for iimgs = 1 : nimgs
    
    [~,name,ext] = fileparts(imlist(iimgs).name);
    %% load images
    img = im2uint8(imread([path name ext]));
%     Data{iimgs}.width_orig  = size(img,2);
%     Data{iimgs}.height_orig = size(img,1);
    %% scale image to smaller size
    img = imresize(img,0.25);
    if options.faceDetected ==0
        [bbox, ~] = facedetect(img);
        save([path name '_bbox.mat'], 'bbox');
    end
    load ([path name '_bbox.mat']);
    bb= step(faceDetector, img);
    if isempty(bb)
        disp('No face detected!');
    else
    %if multiple boxes are detected, choose the largest one
    if size(bb,1) > 1 
        [~,idx] = max(bb(:,3));
        bb = bb(idx,:);
    end
    
    bb = enlargingbbox(bb,1)
    figure(1); imshow(img); hold on;
        rectangle('Position',bb(1,:),'LineWidth',5,'LineStyle','-','EdgeColor','r');
%         draw_shape(Data{iimgs}.shape_gt(:,1),...
%             Data{iimgs}.shape_gt(:,2),'y');
        hold off;
        pause;
    end
     Data{iimgs}.imgCrop = rgb2gray(imcrop(img,bbox));
    Data{iimgs}.width  = ceil(bbox(3));
    Data{iimgs}.height = ceil(bbox(4));
    %% load shape parameter
%     fid = fopen([dbpath_data imlist(iimgs).name(3:end - 4) 'Para.data'],'rb');
%     Data{iimgs}.para_gt = fread(fid,inf,'float');
%     fclose(fid);

%      Data{iimgs}.para_gt = NaN; %if No ground truth
%     %% load transformation parameter
%     %fid = fopen([dbpath_data imlist(iimgs).name(3:end - 4) 'Trans.data'],'rb');
     load([path name '_ground_truth.mat']);
     Data{iimgs}.gtID = wID;
     Data{iimgs}.gtEP = wEP;
     trans = zeros(6,1);
     trans(1:3) = Rot2Ang(R).*(180/pi);
     trans(4:6) = T;
     Data{iimgs}.gtTrans = trans;
     
     load([path name '_campar.mat']);
     
%     trans = zeros(7,1);
%     trans(1:3) = Rot2Ang(Rc).*(180/pi); %transfer rotation matrix to angle degree
%     trans(4:6) = Tc;
%     trans(7) = fx;
%     Data{iimgs}.trans_gt = trans;
    
    Data{iimgs}.centerX = 0.25*cx - bbox(1);
    Data{iimgs}.centerY = 0.25*cy - bbox(2);
    Data{iimgs}.focalX = 0.25*fx;
    Data{iimgs}.focalY = 0.25*fy;
    Data{iimgs}.Rc = Rc;
    Data{iimgs}.Tc = Tc;
    
    %fclose(fid);


    %Data{iimgs}.isdet = 0;
     %Data{iimgs}.img = img;
    
end

end

function region = enlargingbbox(bbox, scale)

region(1) = floor(bbox(1) - (scale - 1)/2*bbox(3));
region(2) = floor(bbox(2) - (scale - 1)/2*bbox(4));

region(3) = floor(scale*bbox(3));
region(4) = floor(scale*bbox(4));

end
