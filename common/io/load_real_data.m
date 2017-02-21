function Data = load_real_data(dbpath_img, dbpath_data, options)
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

slash = options.slash;
dbname = options.datasetName;

imlist = dir([dbpath_img '*.*g']);
nimgs  = length(imlist);
%nimgs  = 10;

%% face detection
% Create a cascade detector object.
%faceDetector = vision.CascadeObjectDetector();

isdetected     = zeros(length(nimgs), 1);

Data = cell(nimgs, 1);

for iimgs = 1 : nimgs
    
    %% load images
    img = im2uint8(imread([dbpath_img imlist(iimgs).name]));
    Data{iimgs}.width_orig  = size(img,2);
    Data{iimgs}.height_orig = size(img,1);
    
    %% load shape parameter
%     fid = fopen([dbpath_data imlist(iimgs).name(3:end - 4) 'Para.data'],'rb');
%     Data{iimgs}.para_gt = fread(fid,inf,'float');
%     fclose(fid);

%      Data{iimgs}.para_gt = NaN; %if No ground truth
%     %% load transformation parameter
%     %fid = fopen([dbpath_data imlist(iimgs).name(3:end - 4) 'Trans.data'],'rb');
     load([dbpath_img imlist(iimgs).name(1:end - 6) 'cam.mat']);
     Data{iimgs}.centerx = camera.cx;
     Data{iimgs}.centery = camera.cy;
     
     load([dbpath_img imlist(iimgs).name(1:end - 6) 'Trans.mat']);
     Data{iimgs}.trans_gt = trans;
     
     if exist([dbpath_img imlist(iimgs).name(1:end - 6) 'Para.mat'])
         load([dbpath_img imlist(iimgs).name(1:end - 6) 'Para.mat']);
         Data{iimgs}.para_gt = para;
     else
         Data{iimgs}.para_gt = NaN;
     end
%     trans = zeros(7,1);
%     trans(1:3) = Rot2Ang(Rc).*(180/pi); %transfer rotation matrix to angle degree
%     trans(4:6) = Tc;
%     trans(7) = fx;
    %fclose(fid);
    %Data{iimgs}.isdet = 0;
     Data{iimgs}.img_gray = rgb2gray(img);
    
end

end