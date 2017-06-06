function [Regression,storageNewInitPara,storageNewInitTran,rms, bv] = learn_single_regressor...
    (model, faces, landmarks, bv, Data, newInitShape, newInitTran,options )

%%%%%%%% properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Using randomly ground-truth of image as initial shape for others.
%% 2. Using multi-scale images.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nData = length(Data);

% if flipping data
if options.flipFlag == 1
    nData = nData * 2;
end

innerLandIdx = landmarks.inner;
outerLandIdxSet= landmarks.boundary;
%% the fixed parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
para_dim  = options.paraSize;%size(ShapeModel.MeanShape,1);
%MeanShape2 = zeros(207);%vec_2_shape(ShapeModel.MeanShape);
n_points   = size(innerLandIdx,2);%+size(outerLandIdxSet,1);
if options.useBoundary ==1 
    n_points   = size(innerLandIdx,2)+size(outerLandIdxSet,1);
end
n_cascades = options.n_cascades;

current_cascade = options.current_cascade;
n_init_randoms  = options.n_init_randoms;
descSize       = options.descSize;
descBins       = options.descBins;

if strcmp(options.descType,'xx_sift') == 1
    desc_dim        = 8 * descBins * descBins; % xx_sift
elseif strcmp(options.descType,'vlsift') == 1
    desc_dim        = 128;
elseif strcmp(options.descType,'hog') == 1
    desc_dim        = 2 * 2 * 31; % hog
else
    desc_dim        = options.descRawWin * options.descRawWin;%raw
end

%% initial matrices used for storing descriptors and delta shape %%%%%%%%%%
storageInitPara = zeros(nData*n_init_randoms, para_dim);
storageGtPara   = zeros(nData*n_init_randoms, para_dim);
storageNewInitTran = cell(nData*n_init_randoms, 1);
storageNewInitPara = zeros(nData*n_init_randoms, para_dim);
storageInitDesc  = zeros(nData*n_init_randoms,desc_dim*n_points);
storageDelPara  = zeros(nData*n_init_randoms, para_dim);
%storage_bbox       = zeros(nData*n_init_randoms,4);

%% set the current canvas size for multi-scale feature
% currentScale = cascade_img_scale(options.scaleFactor,current_cascade,n_cascades);

for idata = 1 : nData       
    %% the information of i-th image
    %disp(Data(idata).img);
    CrTensor = tensor(model.Cr);
    disp(['Stage: ' num2str(options.current_cascade) ' - Image: ' num2str(idata)]);  
    img   = Data{idata}.imgCrop;
    trueID = Data{idata}.gtID;
    trueEP = Data{idata}.gtEP;
    landmark2d = Data{idata}.lm2d';
    lmFlag = Data{idata}.lmFlag;
    
%     %% for test
%     landmarksFlag = ones(1,66);
%     for i = 1:66
%         x = landmark2d(1,i);
%         y = landmark2d(2,i);
%         if (x<0 || x >size(img,1) || y < 0 || y > size(img,2))
%             landmarksFlag(i) = 0;
%         end;
%     end
%     lambdaEP = 100;
%     lambdaID = 100;
%     [R,T,s,wID,wEP] = landmark_fitting(landmark2d,landmarksFlag, model, landmarks,lambdaEP,lambdaID,faces);
%     %% Calculate R,T,s for weak perspective projection
%     meanshape = coeff2shape(CrTensor,mean(model.Wids),mean(model.Weps));
%     %meanshape = reshape(meanshape,3,[]);
%     landmark3d = meanshape(innerLandIdx,:);
%     [R,T,s] = WeakPerspective(landmark2d(1:49,:)',landmark3d');
    
    
    %% ==
%     shape = coeff2shape(CrTensor,wID,wEP);
%     landmark3d = shape(:,innerLandIdx);
%     Proj = cal_weak_perspective(shape, s, R, T);
%     figure(1);imshow(img);hold on;
%     plot(Proj(1,:),Proj(2,:),'r+');
%     plot(landmark2d(1,:),landmark2d(2,:),'g*');hold off;
%     pause;
%     
%     ground = coeff2shape(CrTensor,paraID,paraEP);
%     keyproj = cal_weak_perspective(ground(:,innerLandIdx),s,R,T);
%     figure(2);imshow(img);hold on;
%     proj = cal_weak_perspective(ground,s,R,T);
%     plot(proj(1,:),proj(2,:),'r+');hold on;
%     plot(keyproj(1,:),keyproj(2,:),'b*');hold on;
%     plot(landmark2d(1,:),landmark2d(2,:),'g*');
%     hold off;
%     pause;
%     figure(5);
%     rp = defrp;
%         rp.phi = 0; % frontal face
%         rp.width = 500;
%         rp.height = 500;
%         rp.theta = 0.5*pi;
%         rp.alpha = pi;
% %         mesh_world_coord_3_N = get_face_template_world_coord(R,T,wID,wEP,Bilinear);
% % target_3_N = target_N_3';
% % shape = [mesh_world_coord_3_N(:); target_3_N(:)];
% % tex = [200 * repmat([1;1;1],[1,size(mesh_world_coord_3_N,2)]) 200 * repmat([1;0;0],[1,size(target_3_N,2)])];
% % facesall = [faces;f + size(mesh_world_coord_3_N,2)];
% % display_face(shape, tex, facesall, rp);
%         
%         ground = coeff2shape(CrTensor,paraID,paraEP);
%         est = coeff2shape(CrTensor,wID,wEP);
%         
%         tex = [200 * repmat([1;1;1],[1,size(ground,2)]) 200 * repmat([1;0;1],[1,size(est,2)])];
%         shape = [ground(:);est(:)];
%         facesall = [faces; faces+size(ground,2)];
%      
%         display_face(shape, tex, facesall, rp);
%         pause;
%     
    %% if the first cascade
    if ( current_cascade == 1 )        
%         % randomize which shape is used for initial position
%         rIdx = randi([1,nData],n_init_randoms);
%         % iterations of n initial points
        for ir = 1 : n_init_randoms
            % get random positions and inital shape indexs
%             idx    = rIdx(ir);
%             initID = Data{idx}.gtID; %% get randomly para from other samples        
        %% for test
%             landmarksFlag = ones(1,66);
%             for i = 1:66
%                 x = landmark2d(1,i);
%                 y = landmark2d(2,i);
%                 if (x<0 || x >size(img,1) || y < 0 || y > size(img,2))
%                     landmarksFlag(i) = 0;
%                 end;
%             end
            lambdaEP = 100;
            lambdaID = 100;
            [R,T,s,wID,wEP] = landmark_fitting(landmark2d,lmFlag, model, landmarks,lambdaEP,lambdaID,faces);

            initID = wID;
            initEP = wEP;

            randIndex = randi([1,nData],1);
            initID = Data{randIndex}.gtID;
            initEP = Data{randIndex}.gtEP;
            % initID = mean(model.Wids);
            % initEP = mean(model.Weps);
%% scale coarse to fine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             cropImScale = imresize(img,currentScale);
%             initID = initID * currentScale;
%             trueID    = paraID      * currentScale;      
            %% compute the descriptors and delta_shape %%%%%%%%%%%%%%%%%%%%            
            % storing the initial shape
            storageInitPara((idata-1)*n_init_randoms+ir,:) = [initID initEP];
            %% calculate all key projections
            reconShape = coeff2shape(CrTensor,initID,initEP);
            %% update R,T,s         
            validBoundary = outerLandIdxSet(logical(lmFlag(50:66)),:);
            validLmInner = innerLandIdx(logical(lmFlag(1:49)));  
            target_2_K = landmark2d(:,logical(lmFlag));
            tmp = R * reconShape + repmat(T,[1 size(reconShape,2)]);
            validLmBoundary = get_boundary_vertex(tmp,faces,validBoundary);
            validLm = [validLmInner validLmBoundary];
            source_3_K = reconShape(:,validLm);
            [R,T,s] = weak_perspective(target_2_K,source_3_K);
            clear tmp;
            
            allKey = innerLandIdx;
            if options.useBoundary ==1
                tmp= get_boundary_vertex(reconShape, faces, outerLandIdxSet);
                bv((idata-1)*n_init_randoms+ir,:) = tmp;
                allKey = [innerLandIdx tmp];
            end           
            keyPoints = reconShape(:,allKey);
            initProj = cal_weak_perspective(keyPoints, s, R, T);
            %% show result
            if options.debugMode ==1
                groundShape = coeff2shape(CrTensor,trueID , trueEP);
                keyPoints = groundShape(:,allKey);
                gtProj = cal_weak_perspective(keyPoints, s,R,T);
                figure; imshow(img); 
                hold on;plot(initProj(1,:),initProj(2,:),'r.');
                hold on;plot(gtProj(1,:),gtProj(2,:),'g.');
                hold off;title(['Projection:iter' num2str(current_cascade)]);  
                pause;
            end
            %%  storing the the descriptors
            tmp = local_descriptors(img, initProj, descSize, descBins, options );  
            if options.useBoundary ==1
                tmp(end-17:end,:) = tmp(end-17:end,:)*0.5; %assign weight 0.5 to face boundary
            end
            tmp(1:10,:) = tmp(1:10,:)*0; % do not use eyebrows
            storageInitDesc((idata-1)*n_init_randoms+ir,:) = tmp(:);
            %% storing delta shape
            paraResidual = [initID initEP] - [trueID trueEP];
            storageDelPara((idata-1)*n_init_randoms+ir,:) = paraResidual;
            storageGtPara((idata-1)*n_init_randoms+ir,:) = [trueID trueEP];   
            storageNewInitTran{(idata-1)*n_init_randoms+ir}.R = R;
            storageNewInitTran{(idata-1)*n_init_randoms+ir}.T = T;
            storageNewInitTran{(idata-1)*n_init_randoms+ir}.s = s;
        end
        
    else 
        % for higher cascaded levels
        for ir = 1 : n_init_randoms  
            initPara = newInitShape((idata-1)*n_init_randoms+ir,:);
            initID = initPara(1:50);
            initEP = initPara(51:end);
            R = newInitTran{(idata-1)*n_init_randoms+ir}.R;
            T = newInitTran{(idata-1)*n_init_randoms+ir}.T;
            
            %% scale coarse to fine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             cropImScale = imresize(img,currentScale);
%             initID = initID * currentScale;
%             trueID = trueID      * currentScale;          
            %% compute the descriptors and delta_shape %%%%%%%%%%%%%%%%%%%% 
            % storing the initial shape
            storageInitPara((idata-1)*n_init_randoms+ir,:) = [initID initEP];           
            %% calculate all key projections
            reconShape = coeff2shape(CrTensor,initID,initEP);
            allKey = innerLandIdx;         
            if options.useBoundary ==1
                tmp= get_boundary_vertex(reconShape, faces, outerLandIdxSet);
                bv((idata-1)*n_init_randoms+ir,:) = tmp;
                allKey = [innerLandIdx tmp];
            end         
            keyPoints = reconShape(:,allKey);
            %% update R,T,s         
            validBoundary = outerLandIdxSet(logical(lmFlag(50:66)),:);
            validLmInner = innerLandIdx(logical(lmFlag(1:49)));  
            target_2_K = landmark2d(:,logical(lmFlag));
            tmp = R * reconShape + repmat(T,[1 size(reconShape,2)]);
            validLmBoundary = get_boundary_vertex(tmp,faces,validBoundary);
            validLm = [validLmInner validLmBoundary];
            source_3_K = reconShape(:,validLm);
            [R,T,s] = weak_perspective(target_2_K,source_3_K);
            clear tmp;
            
            initProj = cal_weak_perspective(keyPoints, s, R, T);
            %% show results
            if options.debugMode ==1
                groundShape = coeff2shape(CrTensor,trueID ,trueEP);
                keyPoints = groundShape(:,allKey);
                gtProj = cal_weak_perspective(keyPoints, s,R,T);
                figure; imshow(img); 
                hold on;plot(initProj(1,:),initProj(2,:),'r.');
                hold on;plot(gtProj(1,:),gtProj(2,:),'g.');
                hold off;title(['Projection:iter' num2str(current_cascade)]);     
            end
            
            tmp = local_descriptors( img, initProj, descSize, descBins, options);
            if options.useBoundary ==1
                tmp(end-17:end,:) = tmp(end-17:end,:)*0.5;%assign weight 0.5 to face boundary
            end
            tmp(1:10,:) = tmp(1:10,:)*0;% do not use eyebrows
            storageInitDesc((idata-1)*n_init_randoms+ir,:) = tmp(:);       
            % storing delta shape
            paraResidual = [initID initEP] - [trueID trueEP];
            storageDelPara((idata-1)*n_init_randoms+ir,:) = paraResidual;  
            storageGtPara((idata-1)*n_init_randoms+ir,:) = [trueID trueEP];   
            storageNewInitTran{(idata-1)*n_init_randoms+ir}.R = R;
            storageNewInitTran{(idata-1)*n_init_randoms+ir}.T = T;
            storageNewInitTran{(idata-1)*n_init_randoms+ir}.s = s;
            
        end       
    end
  
    clear img;
    
end

%% solving multivariate linear regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('solving linear regression problem...');
Regression = linreg( storageInitDesc, storageDelPara, ...
    options.lambda(current_cascade) );

%% updading the new shape %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('updadting the para...');
del_shape = regress( storageInitDesc, Regression );
nsamples = size(storageInitDesc,1);

for isample = 1 : nsamples
  
    truePara      = storageInitPara(isample,:) - del_shape(isample,:);
%     trueID      = trueID / currentScale;
    storageNewInitPara(isample,:) = truePara;
    
end

%% compute errors
err = zeros(nsamples,1);

for i = 1:nsamples
    pr_para = storageNewInitPara(i,:);
    gt_para = storageGtPara(i,:);
    err(i) = rms_err( pr_para, gt_para, options);  
end

rms = 100*mean(err);
disp(['ERR average: ' num2str(100*mean(err))]);

clear storage_init_para;
clear storage_gt_para;
clear storage_init_desc;
clear storage_del_para;

