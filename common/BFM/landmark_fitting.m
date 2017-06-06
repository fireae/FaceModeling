function [R,T,s,wID,wEP] = landmark_fitting(landmark2d,landmarksFlag, model, landmarks3d,lambdaEP,lambdaID,faces)
    
    % Constants
    scale2mm = 100;   
    validBoundary = landmarks3d.boundary(logical(landmarksFlag(50:66)),:);
    
    validLmInner = landmarks3d.inner(logical(landmarksFlag(1:49)));
    % 24 landmarks
    %target_K_2 = landmark2d(logical(landmarksFlag),:);
    target_2_K = landmark2d(:,logical(landmarksFlag));
    
    MaxItr = 2;
    CrTensor = tensor(model.Cr);   % reduced core tensor
    err = 4;
    wID = mean(model.Wids);
    wEP = mean(model.Weps);
    R = eye(3); T = zeros(3,1);
    for itr=1:MaxItr
        %disp(['Iteration: ' num2str(itr) '/' num2str(MaxItr)]);
        shape = coeff2shape(CrTensor,wID,wEP);
        shape = R * shape + repmat(T,[1 size(shape,2)]);
        validLmBoundary = get_boundary_vertex(shape,faces,validBoundary);
        validLm = [validLmInner validLmBoundary];
        % R and T
        source_3N_1 = scale2mm * double(ttm(CrTensor,{wID,wEP},[2,3]));
        source_3_N = reshape(source_3N_1,3,[]);
        source_3_K = source_3_N(:,validLm);
        [R,T,s] = weak_perspective(target_2_K,source_3_K);
        
        % Expression
        wEP = cal_expression(R,T,s,wID,target_2_K,model,validLm,lambdaEP);
%         wID = mean(model.Wids);
%         wEP = mean(model.Weps);
%         figure(3);
%         rp = defrp;
%         rp.phi = 0; % frontal face
%         rp.width = 500;
%         rp.height = 500;
%         rp.theta = 0.5*pi;
%         rp.alpha = pi;
%         reconShape = coeff2shape(CrTensor,wID,wEP);
%         
%         
%         tex = 200 * repmat([1;1;1],[1,size(reconShape,2)]);
%         shape = reconShape(:);
%         display_face(shape, tex, faces, rp);
%         pause;
        % Identity
        wID = cal_identity(R,T,s,wEP,target_2_K,model,validLm,lambdaID);
%         figure(4);
%         rp = defrp;
%         rp.phi = 0; % frontal face
%         rp.width = 500;
%         rp.height = 500;
%         rp.theta = 0.5*pi;
%         rp.alpha = pi;
%         reconShape = coeff2shape(CrTensor,wID,wEP);        
%         tex = 200 * repmat([1;1;1],[1,size(reconShape,2)]);
%         shape = reconShape(:);
%         display_face(shape, tex, faces, rp);
%         pause;
        %err = show_fitting_error(R,T,wID,wEP,model,faces,landmark2d,f, err);
    end

end

function  wID = cal_identity(R,T,s,wEP,validLandmarks2d,model,validLandmarksModel,lambdaID)

    % Constants
    scale2mm = 100;
    wIDmean = mean(model.Wids);
    wIDstd  = std(model.Wids);
    
    CrTensor = tensor(model.Cr);
    Ptmp = scale2mm * double(ttm(CrTensor,{wEP},[3]));
    P_3_N_50 = reshape(Ptmp,3,[],50);
    
    K = length(validLandmarks2d); % 49 landmarks
    
    % P 
%     idx_3_K = zeros(3,K);
%     idx_3_K(1,:) = 3*(validLandmarksModel-1)+1;
%     idx_3_K(2,:) = 3*(validLandmarksModel-1)+2;
%     idx_3_K(3,:) = 3*(validLandmarksModel-1)+3;
%     idx_3K_1 = idx_3_K(:);
%     P_3K_50 = P_3N_50(idx_3K_1,:);
    
    P_3_K_50 = P_3_N_50(:,validLandmarksModel,:);
    RP_2_K_50 = zeros(2,K,50);
    for i = 1:25
        RP_2_K_50(:,:,i) = R(1:2,:)*P_3_K_50(:,:,i);
    end
    sRP_2_K_50 = s.*RP_2_K_50;
    sRP_2K_50 = reshape(sRP_2_K_50,[],50);
    
    % (T-D)    
    t = T(1:2);
    D_2_K = validLandmarks2d;
    t_2_K = repmat(t,[1,K]);
    t_sub_D = t_2_K - D_2_K;
    t_sub_D_2K_1 = reshape(t_sub_D,[],1);
    
    %Q
    id_std = wIDstd;
    id_var = id_std.*id_std;
    Qa = eye(50);
    Qa(1:50+1:end) = 1./id_var(1:50);
    
    % Ax = b
    % A = P'*P+lambda*Q
    % b = lambda*Q*mean - P'*(T-D)
    A = sRP_2K_50'*sRP_2K_50 + lambdaID * Qa;
    b = lambdaID * Qa * wIDmean' - sRP_2K_50' * t_sub_D_2K_1;
    x = A\b;
    wID = x';
end

function  wEP = cal_expression(R,T,s,wID,validLandmarks2d,model,validLandmarksModel,lambdaEP)

    % Constants
    scale2mm = 100;
    wEPmean = mean(model.Weps);
    wEPstd  = std(model.Weps);
    
    CrTensor = tensor(model.Cr);
    Ptmp = scale2mm * double(ttm(CrTensor,{wID},[2]));
    P_3_N_25 = reshape(Ptmp,3,[],25);
    K = length(validLandmarks2d); % 24 landmarks
    
    % P:s.*R*PCA
%     idx_3_K = zeros(3,K);
%     idx_3_K(1,:) = 3*(validLandmarksModel-1)+1;
%     idx_3_K(2,:) = 3*(validLandmarksModel-1)+2;
%     idx_3_K(3,:) = 3*(validLandmarksModel-1)+3;
%     idx_3K_1 = idx_3_K(:);
    P_3_K_25 = P_3_N_25(:,validLandmarksModel,:);
    RP_2_K_25 = zeros(2,K,25);
    for i = 1:25
        RP_2_K_25(:,:,i) = R(1:2,:)*P_3_K_25(:,:,i);
    end
    sRP_2_K_25 = s.*RP_2_K_25;
    sRP_2K_25 = reshape(sRP_2_K_25,[],25);
    
    % T_sub_Li: t - landmark2d
    t = T(1:2);
    D_2_K = validLandmarks2d;
    t_2_K = repmat(t,[1,K]);
    t_sub_D = t_2_K - D_2_K;
    t_sub_D_2K_1 = reshape(t_sub_D,[],1);
    
    %Q
    ep_std = wEPstd;
    ep_var = ep_std.*ep_std;
    Qb = eye(25);
    Qb(1:25+1:end) = 1./ep_var(1:25);
    
    % Ax = b
    % A = P'*P+lambda*Q
    % b = lambda*Q*mean - P'*R'*(T-D)
    A = sRP_2K_25'*sRP_2K_25 + lambdaEP * Qb;
    b = lambdaEP * Qb * wEPmean' - sRP_2K_25' * t_sub_D_2K_1;
    x = A\b;
    wEP = x';
end

