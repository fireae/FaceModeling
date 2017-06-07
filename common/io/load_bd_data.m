function Data = load_bd_data( path, options )
%LOAD_BD_DATA Summary of this function goes here
%   Detailed explanation goes here
facelist = dir([path '*data']);
nfaces  = length(facelist);
%nimgs  = 10;

Data = cell(nfaces, 1);
Label = cell(nfaces,1);
validFace = zeros(1,nfaces);
for ifaces = 1 : nfaces
    
    [~,name,ext] = fileparts(facelist(ifaces).name);
    [tmp1 ctg sub1 sub2] = strread(name, '%s %s %s %s', 'delimiter', '_');
        if (strcmp(char(ctg(1)),'O')) 
         %||(strcmp(char(ctg(1)),'CR')) ...
         %||(strcmp(char(ctg(1)),'PR')) ...
         %(strcmp(char(ctg(1)),'YR')&strcmp(char(sub1(1)),'L45')) ...
         %|| (strcmp(char(ctg(1)),'YR')&strcmp(char(sub1(1)),'R45')) ...
         %(strcmp(char(ctg(1)),'YR')&strcmp(char(sub1(1)),'L90')) ...
         %||(strcmp(char(ctg(1)),'YR')&strcmp(char(sub1(1)),'R90')) ...
         %||(strcmp(char(ctg(1)),'YR')&strcmp(char(sub1(1)),'L90')) ...
         %||(strcmp(char(ctg(1)),'YR')&strcmp(char(sub1(1)),'R90')) ...
         validFace(ifaces) = 0;
        Data{ifaces}.flag =0;
        else 
        Data{ifaces}.flag = 1;
        validFace(ifaces) =1;

    Label{ifaces}.name = name;
    %%load 2D landmarks
    fid = fopen([path name ext],'rb');
    landmarks = fread(fid,inf,'float');
    fclose(fid);
    Data{ifaces}.lm2d = reshape(landmarks,2,[])';
    
    %% load images
    img = im2uint8(imread([path name '.png']));

    Data{ifaces}.width_orig = size(img,2);
    Data{ifaces}.height_orig = size(img,1);
    %% crop image to normal size
    bbox = getbbox(Data{ifaces}.lm2d);
    
    %% enlarge region of face
    region     = enlargingbbox(bbox, 1.2);
    region(2)  = double(max(region(2), 1));
    region(1)  = double(max(region(1), 1));
    
    bottom_y   = double(min(region(2) + region(4) - 1, ...
        Data{ifaces}.height_orig));
    right_x    = double(min(region(1) + region(3) - 1, ...
        Data{ifaces}.width_orig));
    
    imgCrop = img(region(2):bottom_y, region(1):right_x, :);
    %% normalize image
            %tmp = reshape(double(imgCrop),1,[]);
            %Imax = max(tmp);
            %Imin = min(tmp);
            %imgCrop = reshape(uint8(255*(tmp - Imin)/(Imax - Imin)),size(imgCrop));    
    %% recalculate the location of groundtruth shape and bounding box
    Data{ifaces}.lm2d = bsxfun(@minus, Data{ifaces}.lm2d,...
        double([region(1) region(2)]));
    Data{ifaces}.imgCrop = rgb2gray(imgCrop);
    Data{ifaces}.width  = size(img,2);
    Data{ifaces}.height = size(img,1);
    
    landmarksFlag = ones(1,66);
            for i = 1:66
                x = Data{ifaces}.lm2d(i,1);
                y = Data{ifaces}.lm2d(i,2);
                if (x<0 || x >size(imgCrop,1) || y < 0 || y > size(imgCrop,2))
                    landmarksFlag(i) = 0;
                end;
            end
    Data{ifaces}.lmFlag = landmarksFlag;
    %% load shape parameter
     load([path name '_ground_truth.mat']);
     Data{ifaces}.gtID = wID;
     Data{ifaces}.gtEP = wEP;
     
     if options.debugMode ==1
         figure(1);
         imshow(imgCrop);
         hold on; 
         plot(Data{ifaces}.lm2d(:,1),Data{ifaces}.lm2d(:,2),'g*');
         hold off;
         pause;
     end
        end
    
        
end
Data = Data(logical(validFace));
Label = Lable(logical(validFace));
save([options.ResultPath 'Lable.mat','Label');
end

function region = enlargingbbox(bbox, scale)

region(1) = floor(bbox(1) - (scale - 1)/2*bbox(3));
region(2) = floor(bbox(2) - (scale - 1)/2*bbox(4));

region(3) = floor(scale*bbox(3));
region(4) = floor(scale*bbox(4));

end

