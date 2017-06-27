function [Regression,storageNewInitPara,storageNewInitTran,rms] = learn_single_regressor...
    (model, faces, landmarks, Data, newInitShape, newInitTran,options )

nData = length(Data);
% if flipping data
if options.flipFlag == 1
    nData = nData * 2;
end

innerLandIdx = landmarks.inner;
outerLandIdxSet= landmarks.boundary;
%% the fixed parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
para_dim  = options.paraSize;%size(ShapeModel.MeanShape,1);

n_points   = size(innerLandIdx,2);%+size(outerLandIdxSet,1);
if options.useBoundary ==1 
    n_points   = size(innerLandIdx,2)+size(outerLandIdxSet,1);
end

current_cascade = options.current_cascade;
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
storageInitPara = zeros(nData, para_dim);
storageGtPara   = zeros(nData, para_dim);
storageNewInitTran = cell(nData, 1);
storageNewInitPara = zeros(nData, para_dim);
storageInitDesc  = zeros(nData,(2+desc_dim)*n_points);
storageDelPara  = zeros(nData, para_dim);
CrTensor = tensor(model.Cr);
if ( current_cascade == 1 )    
 parfor idata = 1 : nData       
    %% the information of i-th image
    %disp(Data(idata).img);
    disp(['Stage: ' num2str(current_cascade) ' - Image: ' num2str(idata)]);  
    img   = Data{idata}.imgCrop;
    trueID = Data{idata}.gtID;
    trueEP = Data{idata}.gtEP;
    landmark2d = Data{idata}.lm2d';
    lmFlag = Data{idata}.lmFlag;
      
            lambdaEP = 100;
            lambdaID = 100;
            [R,T,s,wID,wEP] = landmark_fitting(landmark2d,lmFlag, model, landmarks,lambdaEP,lambdaID,faces);

            initID = wID;
            initEP = wEP;
            initID = mean(model.Wids);
            initEP = mean(model.Weps); 
        %% compute the descriptors and delta_shape %%%%%%%%%%%%%%%%%%%%            
            % storing the initial shape
            storageInitPara(idata,:) = [initID initEP];
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
                        
            allKey = innerLandIdx;
            if options.useBoundary ==1
                tmp= get_boundary_vertex(reconShape, faces, outerLandIdxSet);
                bv(idata,:) = tmp;
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
            %% uv between landmarks and projections
            uv = initProj - landmark2d;
            %%
            storageInitDesc(idata,:) = [reshape(tmp,1,[]) reshape(uv,1,[])];
            %storageInitDesc(idata,:) = tmp(:);
            %% storing delta shape
            paraResidual = [initID initEP] - [trueID trueEP];
            storageDelPara(idata,:) = paraResidual;
            storageGtPara(idata,:) = [trueID trueEP];   
            storageNewInitTran{idata}.R = R;
            storageNewInitTran{idata}.T = T;
            storageNewInitTran{idata}.s = s;
 end
        
    else 
        % for higher cascaded levels
        parfor idata = 1 : nData    
            disp(['Stage: ' num2str(current_cascade) ' - Image: ' num2str(idata)]);  
            img   = Data{idata}.imgCrop;
            trueID = Data{idata}.gtID;
            trueEP = Data{idata}.gtEP;
            landmark2d = Data{idata}.lm2d';
            lmFlag = Data{idata}.lmFlag;
    
            initPara = newInitShape(idata,:);
            initID = initPara(1:50);
            initEP = initPara(51:end);
            R = newInitTran{idata}.R;
            T = newInitTran{idata}.T;
            
         %% compute the descriptors and delta_shape %%%%%%%%%%%%%%%%%%%% 
            % storing the initial shape
            storageInitPara(idata,:) = [initID initEP];           
            %% calculate all key projections
            reconShape = coeff2shape(CrTensor,initID,initEP);
            allKey = innerLandIdx;         
            if options.useBoundary ==1
                tmp= get_boundary_vertex(reconShape, faces, outerLandIdxSet);
                bv(idata,:) = tmp;
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
            %% uv between landmarks and projections
            uv = initProj - landmark2d;
            %%
            storageInitDesc(idata,:) = [reshape(tmp,1,[]) reshape(uv,1,[])];
            %storageInitDesc(idata,:) = tmp(:);       
            % storing delta shape
            paraResidual = [initID initEP] - [trueID trueEP];
            storageDelPara(idata,:) = paraResidual;  
            storageGtPara(idata,:) = [trueID trueEP];   
            storageNewInitTran{idata}.R = R;
            storageNewInitTran{idata}.T = T;
            storageNewInitTran{idata}.s = s;
        end
        
end
  
% end

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
    storageNewInitPara(isample,:) = truePara;
    
end

%% compute errors
err = zeros(nsamples,1);

for i = 1:nsamples
    pr_para = storageNewInitPara(i,:);
    gt_para = storageGtPara(i,:);
    err(i) = rms_err( pr_para, gt_para, options);  
end

rms = mean(err);
disp(['ERR average: ' num2str(mean(err))]);

clear storage_init_para;
clear storage_gt_para;
clear storage_init_desc;
clear storage_del_para;

