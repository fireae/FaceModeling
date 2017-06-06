function [aligned_shape,R,T,s] = cascaded_regress (landmark2d,lmFlag,CrTensor, faces, Landmarks,...
    LearnedCascadedModel, img, initID, initEP, index,R,T,s,options )


%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nCascades = LearnedCascadedModel{1}.n_cascades;
descSize  = LearnedCascadedModel{1}.descSize;
descBins  = LearnedCascadedModel{1}.descBins;
factor     = options.scaleFactor;

innerLandIdx = Landmarks.inner;
outerLandIdxSet= Landmarks.boundary;

%% iterations of cascades %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ic = 1 : nCascades
    
    current_scale = cascade_img_scale(factor, ic, nCascades);
    
    options.current_cascade = ic;
    
    cropImScale = imresize(img,current_scale);
    initID   = initID * current_scale;
    
    reconShape = coeff2shape(CrTensor, initID, initEP);
    %% update R,T,s         
            validBoundary = Landmarks.boundary(logical(lmFlag(50:66)),:);
            validLmInner = Landmarks.inner(logical(lmFlag(1:49)));  
            target_2_K = landmark2d(:,logical(lmFlag));
            tmp = R * reconShape + repmat(T,[1 size(reconShape,2)]);
            validLmBoundary = get_boundary_vertex(tmp,faces,validBoundary);
            validLm = [validLmInner validLmBoundary];
            source_3_K = reconShape(:,validLm);
            [R,T,s] = weak_perspective(target_2_K,source_3_K);
            clear tmp;
    
    allKey = innerLandIdx;
    if options.useBoundary ==1
%         if ic == 1          
        bv = get_boundary_vertex(reconShape, faces, outerLandIdxSet);
%         end
        allKey = [innerLandIdx bv];
    end

    keyPoints = reconShape(:,allKey);
    initProj = cal_weak_perspective(keyPoints, s,R,T);
    % extract local descriptors
%     init_projection = CalKeyProj(BFMmodel,init_shape,keypoints);
    if options.debugMode ==1
        figure(2); imshow(cropImScale); hold on;
        plot(initProj(1,:), initProj(2,:),'g.');
        title(['Keypoints:iter' num2str(ic)]);      
        hold off;
        
        rp = defrp;
        rp.phi = 0; % frontal face
        rp.width = 500;
        rp.height = 500;
        rp.theta = 0.5*pi;
        rp.alpha = pi;
        figure(3);
        tex = 200 * repmat([1;1;1],[1,size(reconShape,2)]);
        shape = reconShape(:);
        display_face(shape, tex, faces, rp);
        title('3D face:estimation');
        pause;
        %saveas(gcf,[options.ResultPath num2str(idata) '_shape.jpg']);
    end
%     rp = defrp;
%     rp.phi = 0; % frontal face
%     rp.theta = 0.5*pi;
%     rp.alpha = pi;
%     shape = reconShape(:);
%     tex = 200 * repmat([1;1;1],[1,size(reconShape,2)]);
%     figure(5);
%     display_face(shape, tex, faces, rp);  
%     title('3D face:estimated');
%     saveas(gcf,[options.ResultDataPath num2str(idata) 'iter' num2str(ic) 'estimate3D.jpg']);
    desc = local_descriptors(cropImScale, ...
    initProj, descSize, descBins, options);
    if options.useBoundary ==1 
      desc(end-17:end,:) = desc(end-17:end,:)*0.5;
    end
    desc(1:10,:) = desc(1:10,:)*0;
    % regressing
    %delPara = regress( desc(:)', LearnedCascadedModel{ic}.Regression );
    uv = initProj - landmark2d;
    % regressing
    %delPara = regress( desc(:)', LearnedCascadedModel{ic}.Regression );
    delPara = regress( [reshape(desc(:),1,[]) reshape(uv,1,[])], LearnedCascadedModel{ic}.Regression );
 
    %origin_del = bsxfun(@times, vec_2_shape(del_shape'), bbox(3:4));
    
    % estimate the new shape
    tmp_para = [initID initEP] - delPara;
    
    tmp_para = tmp_para / current_scale;
    
    initID    = tmp_para(1:50);
    initEP   =tmp_para(51:end);
    
    aligned_shape = tmp_para;
    
    
    
end


end

