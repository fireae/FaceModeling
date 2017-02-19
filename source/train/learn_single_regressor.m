function [R,storage_new_init_para,rms, bv] = learn_single_regressor...
    (BFMmodel, innerKeypointIndices, boundIdxSet, bv, Data, new_init_shape, options )

%%%%%%%% properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Using randomly ground-truth of image as initial shape for others.
%% 2. Using multi-scale images.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nData = length(Data);

% if flipping data
if options.flipFlag == 1
    nData = nData * 2;
end

%% the fixed parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
para_dim  = options.paraSize;%size(ShapeModel.MeanShape,1);
%MeanShape2 = zeros(207);%vec_2_shape(ShapeModel.MeanShape);
n_points   = size(innerKeypointIndices,2)+size(boundIdxSet,2);
n_cascades = options.n_cascades;

current_cascade = options.current_cascade;
n_init_randoms  = options.n_init_randoms;
desc_size       = options.descSize;
desc_bins       = options.descBins;

if strcmp(options.descType,'xx_sift') == 1
    desc_dim        = 8 * desc_bins * desc_bins; % xx_sift
elseif strcmp(options.descType,'hog') == 1
    desc_dim        = 2 * 2 * 31; % hog
else
    desc_dim        = options.descRawWin * options.descRawWin;%raw
end

%% initial matrices used for storing descriptors and delta shape %%%%%%%%%%
storage_init_para = zeros(nData*n_init_randoms, para_dim);
storage_gt_para   = zeros(nData*n_init_randoms, para_dim);
storage_new_init_para = zeros(nData*n_init_randoms, para_dim);
storage_init_desc  = zeros(nData*n_init_randoms,desc_dim*n_points);
storage_del_para  = zeros(nData*n_init_randoms, para_dim);
%storage_bbox       = zeros(nData*n_init_randoms,4);

%% set the current canvas size for multi-scale feature
current_scale = cascade_img_scale(options.scaleFactor,current_cascade,n_cascades);
%%load data file path
imgDir = options.trainingImageDataPath;
slash = options.slash;
plist = dir([imgDir slash 'im*.*g']);


for idata = 1 : nData       
    %% the information of i-th image
    %disp(Data(idata).img);
    disp(['Stage: ' num2str(options.current_cascade) ' - Image: ' num2str(idata)]);  
    img   = Data{idata}.img_gray;
    para = Data{idata}.para_gt;
    trans = Data{idata}.trans_gt;
     
    imsize = size(img);
    %% if the first cascade
    if ( current_cascade == 1 )        
        % randomize which shape is used for initial position
        rIdx = randi([1,nData],n_init_randoms);
        % iterations of n initial points
        for ir = 1 : n_init_randoms
            
            % get random positions and inital shape indexs
            idx    = rIdx(ir);
            init_para = Data{idx}.para_gt; %% get randomly para from other samples
            %init_para(end-7:end) = para(end-7:end);%use ground truth transformation
 
            init_trans = trans;
            %% scale coarse to fine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cropIm_scale = imresize(img,current_scale);
            init_para = init_para * current_scale;
            true_para    = para      * current_scale;
            
            %% compute the descriptors and delta_shape %%%%%%%%%%%%%%%%%%%%            
            % storing the initial shape
            storage_init_para((idata-1)*n_init_randoms+ir,:) = init_para;
            all_key = innerKeypointIndices;
            
            reconShape = Reconstruct_face(BFMmodel,init_para,init_trans);
            if options.useBoundary ==1
                tmp= get_boundary_vertex(reconShape, BFMmodel.tl, boundIdxSet);
                bv((idata-1)*n_init_randoms+ir,:) = tmp;
                all_key = [innerKeypointIndices tmp];
            end
            
            initProjection = CalKeyProj(reconShape, init_trans, all_key,imsize);
            
            if 0
                groundShape = Reconstruct_face(BFMmodel, true_para, trans);
                gtProjection = CalKeyProj(groundShape, trans, all_key, imsize);
                figure; imshow(cropIm_scale); 
                hold on;
                plot(initProjection(:,1),initProjection(:,2),'r.');
                hold on;
                plot(gtProjection(:,1),gtProjection(:,2),'g.');
                hold off;
                title(['Projection:iter' num2str(current_cascade)]);
                pause;
            end
            % storing the the descriptors
            tmp = local_descriptors( cropIm_scale, initProjection, desc_size, desc_bins, options );  
            tmp(end-17:end,:) = tmp(end-17:end,:)*0.5;
            tmp(1:10,:) = tmp(1:10,:)*0;
            storage_init_desc((idata-1)*n_init_randoms+ir,:) = tmp(:);
            
            % storing delta shape
            para_residual = init_para - true_para;
            storage_del_para((idata-1)*n_init_randoms+ir,:) = para_residual;
            storage_gt_para((idata-1)*n_init_randoms+ir,:) = para;
                              
        end
        
    else 
        % for higher cascaded levels
        for ir = 1 : n_init_randoms
     
            init_para = new_init_shape((idata-1)*n_init_randoms+ir,:)';
            %% scale coarse to fine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cropIm_scale = imresize(img,current_scale);
            init_para = init_para * current_scale;
            true_para = para      * current_scale;
            
            
            %% compute the descriptors and delta_shape %%%%%%%%%%%%%%%%%%%% 
            % storing the initial shape
            storage_init_para((idata-1)*n_init_randoms+ir,:) = init_para;           
            % storing the the descriptors
            all_key = innerKeypointIndices;         
            reconShape = Reconstruct_face(BFMmodel, init_para, trans);
            if options.useBoundary ==1
                all_key = [innerKeypointIndices bv((idata-1)*n_init_randoms+ir,:)];
            end           
            initProjection = CalKeyProj(reconShape, trans, all_key, imsize);
            
            if 0
                groundShape = Reconstruct_face(BFMmodel,true_para, trans);
                gtProjection = CalKeyProj(groundShape, trans, all_key, imsize);
                figure; imshow(cropIm_scale); 
                hold on;
                plot(initProjection(:,1),initProjection(:,2),'r.');
                hold on;
                plot(gtProjection(:,1),gtProjection(:,2),'g.');
                hold off;
                title(['Projection:iter' num2str(current_cascade)]);
                pause;
            end
            
            tmp = local_descriptors( cropIm_scale, initProjection, desc_size, desc_bins, options);
            tmp(end-17:end,:) = tmp(end-17:end,:)*0.5;
            tmp(1:10,:) = tmp(1:10,:)*0;
            storage_init_desc((idata-1)*n_init_randoms+ir,:) = tmp(:);
            
            % storing delta shape
            para_residual = init_para - true_para;
            storage_del_para((idata-1)*n_init_randoms+ir,:) = para_residual;  
            storage_gt_para((idata-1)*n_init_randoms+ir,:) = para;
            
            
        end       
    end
  
    clear img;
    clear cropIm_scale;
    clear shape;
    
end

%% solving multivariate linear regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('solving linear regression problem...');
R = linreg( storage_init_desc, storage_del_para, ...
    options.lambda(current_cascade) );

%% updading the new shape %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('updadting the para...');
del_shape = regress( storage_init_desc, R );
nsamples = size(storage_init_desc,1);

for isample = 1 : nsamples
  
    para      = storage_init_para(isample,:) - del_shape(isample,:);
    para      = para / current_scale;
    storage_new_init_para(isample,:) = para;
    
end

%% compute errors

err = zeros(nsamples,1);

for i = 1:nsamples
    pr_para = storage_new_init_para(i,:);
    gt_para = storage_gt_para(i,:);
    err(i) = rms_err( pr_para, gt_para, options);  
end

rms = 100*mean(err);
disp(['ERR average: ' num2str(100*mean(err))]);

clear storage_init_para;
clear storage_gt_para;
clear storage_init_desc;
clear storage_del_para;
clear storage_transM;

