function [estID,estEP,R,T,s] = face_alignment( model, faces, landmarks,...
    LearnedCascadedModel, Data, index, options )



%% setup the fixed parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nData = length(trainData);
n_init_randoms_test = options.n_init_randoms_test;
% paraDim  = options.paraSize;

% aligned_shape = zeros(n_init_randoms_test, paraDim);

%% randomize which shape is used for initial position
% rIdx = randi([1,nData],n_init_randoms_test);

%% iterations of n initial points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ir = 1 : n_init_randoms_test
    
    %% get random positions and inital shape indexs
%     idx    = rIdx(ir);
%     initID = trainData{idx}.gtID; %% get randomly shape from others
    img = Data.imgCrop;
    %Trans = Data.gtTrans;
%     EP = Data.gtEP;
    landmark2d = Data.lm2d';
    lmFlag = Data.lmFlag;
    
    %% Get initial R T s wID wEP
    
    lambdaEP = 100;
    lambdaID = 100;
    [R,T,s,wID,wEP] = landmark_fitting(landmark2d,lmFlag, model, landmarks,lambdaEP,lambdaID,faces);
    CrTensor = tensor(model.Cr);
    %% test
%     tmp.ID = wID;
%     tmp.EP = wEP;
%     tmp.R = R;
%     tmp.T = T;
%     tmp.s = s;
%     
%     [error1,error2] = visualize_result(Data,tmp,CrTensor, faces, landmarks, Segmentation, options);
    %initTrans = Data.gtTrans;   
    if options.debugMode ==1
        
        estShape = coeff2shape(CrTensor, wID, wEP); 
        transformedestShape = R*estShape+repmat(T,[1,size(estShape,2)]);
        tmp= get_boundary_vertex(transformedestShape, faces, landmarks.boundary);
        allKey = [landmarks.inner tmp];
        keyPoints = estShape(:,allKey); 
        estProj = cal_weak_perspective(keyPoints, s, R, T);
        figure(2);
        imshow(img);
        hold on;plot(estProj(1, :),estProj(2,:),'r*');
%     hold on; plot(gtProj(1,:),gtProj(2,:),'g*');
        hold on;plot(landmark2d(1,:),landmark2d(2,:),'b*');
        hold off;
        title('Initialization');
    
        rp = defrp;
        rp.phi = 0; % frontal face
        rp.width = 500;
        rp.height = 500;
        rp.theta = 0.5*pi;
        rp.alpha = pi;
        figure(3);
        tex = 200 * repmat([1;1;1],[1,size(estShape,2)]);
        shape = estShape(:);
        display_face(shape, tex, faces, rp);
        title('3D face:estimation');
        pause;
    end
    %% detect landmarks using cascaded regression
    
    [estPara, R, T,s] = cascaded_regress( landmark2d,lmFlag,CrTensor, faces, landmarks, ...
        LearnedCascadedModel, img, wID, wEP,index,R, T,s, options);    
end

if n_init_randoms_test == 1
    estID = estPara(1:50);
    estEP = estPara(51:end);
    
else
    estID = estPara(1:50);
    estEP = estPara(51:end);
end

end
