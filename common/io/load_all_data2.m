function [Data] = load_all_data2 (dbpath_img, dbpath_data, options )

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

imlist = dir([dbpath_img slash '*.*g']);
nimgs  = length(imlist);
%nimgs  = 10;

%% face detection
% Create a cascade detector object.
%faceDetector = vision.CascadeObjectDetector();

isdetected     = zeros(length(nimgs), 1);

Data = cell(nimgs, 1);

for iimgs = 1 : nimgs
    
    %% load images
    img = im2uint8(imread([dbpath_img slash imlist(iimgs).name]));
    Data{iimgs}.width_orig  = size(img,2);
    Data{iimgs}.height_orig = size(img,1);
    
    %% load shape parameter
    fid = fopen([dbpath_data slash imlist(iimgs).name(3:end - 4) 'Para.data'],'rb');
    Data{iimgs}.para_gt = fread(fid,inf,'float');
    fclose(fid);

    %% load transformation parameter
    fid = fopen([dbpath_data slash imlist(iimgs).name(3:end - 4) 'Trans.data'],'rb');
    Data{iimgs}.trans_gt = fread(fid,inf,'float');
    fclose(fid);


    load ([dbpath_data 'cam' imlist(iimgs).name(3:end - 4) '.mat'])
    Data{iimgs}.centerx = camera.cx;
    Data{iimgs}.centery = camera.cy;

    %Data{iimgs}.isdet = 0;
     %Data{iimgs}.img = img;
     Data{iimgs}.img_gray = rgb2gray(img);
    
end


end


function region = enlargingbbox(bbox, scale)

region(1) = floor(bbox(1) - (scale - 1)/2*bbox(3));
region(2) = floor(bbox(2) - (scale - 1)/2*bbox(4));

region(3) = floor(scale*bbox(3));
region(4) = floor(scale*bbox(4));

end

