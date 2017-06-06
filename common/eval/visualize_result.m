function [projdist, rmseParaCenter, rmseParaWhole,rmseCenter,rmseWhole,euDist]...
    =  visualize_result(Data, estimation, CrTensor, faces, Landmarks, Segmentation, options)
    img          = Data.imgCrop; % using grayscale   
%     cx             = Data.centerX;          %center along horizontal direction
%     cy             = Data.centerY;          %center clong vertical direction
%     fx             = Data.focalX;
%     fy             = Data.focalY;
%     Rc            = Data.Rc;
%     Tc            = Data.Tc;
%  

    idxValid = logical(Segmentation.segValid);    
    idxCenter = Segmentation.segCenter;
%% load estimation
    EP = estimation.EP;
    ID = estimation.ID;
    s = estimation.s;
    R = estimation.R;
    T = estimation.T;
    
    %% reconstruct shape with ground truth parameters
    gtShape = coeff2shape(CrTensor, Data.gtID, Data.gtEP);
    allKey = Landmarks.inner; 
    transformedGtShape = R*gtShape+repmat(T,[1,size(gtShape,2)]);
    if options.useBoundary ==1
                tmp= get_boundary_vertex(transformedGtShape, faces, Landmarks.boundary);
                allKey = [Landmarks.inner tmp];
    end   
    tmp= get_boundary_vertex(transformedGtShape, faces, Landmarks.boundary);
    allKey = [Landmarks.inner tmp];
    keyPoints = gtShape(:,allKey);
    gtProj = cal_weak_perspective(keyPoints, s, R, T);
    %% reconstruct shape with estimated parameters
    estShape = coeff2shape(CrTensor, ID, EP); 
    keyPoints = estShape(:,allKey); 
    estProj = cal_weak_perspective(keyPoints, s, R, T);
    %% Normalized Project distance
    projdist = norm_Edist(gtProj, estProj, options);
    %% show image with projection
    figure(1);
    imshow(img);
    hold on;plot(estProj(1, :),estProj(2,:),'r*');
    hold on; %plot(gtProj(1,:),gtProj(2,:),'g*');
    hold on;plot(Data.lm2d(:,1),Data.lm2d(:,2),'b*');
    %show shape mask
%     tmp = cal_key_projection(estShape, cx, cy, fx, fy, Rc, Tc);
%     plot(tmp(:,1),tmp(:,2),'b.',);
    title(['Key points projection - projdist =' num2str(projdist) ]);
    legend('estimate','ground truth');
    hold off;
    %% RMSE with ground truth point cloud
%     vec = Data.pCloud;
%     [Ricp Ticp ER t] = icp(estShape(:,Segmentation.segCenter),vec, 15);
%     rmse2 = ER(end);
%     disp(['Distance with Point Cloud:' num2str(rmse2)]);
%        rmse2 = 0;
%     Dicp = Ricp * estShape + repmat(Ticp, 1, length(IdxValid));
    
    
    %% RMSE
    error = abs(estShape(3,:)-gtShape(3,:));
    rmseParaCenter = sqrt(sum(error(idxCenter).^2)/length(error(idxCenter)));
    rmseParaWhole = sqrt(sum(error(idxValid).^2)/length(error(idxValid)));
    disp(['rmse with gt parameters: Face:' num2str(rmseParaCenter) 'Head:' num2str(rmseParaWhole)]);
    %% show 3D error map
    %% RMSE with rigistered mesh
    mesh = Data.mesh;
    estShape = Data.worldR*estShape + repmat(Data.worldT,[1 size(estShape,2)]);
    error = abs(estShape(3,:)-mesh(3,:));
    rmseCenter = sqrt(sum(error(idxCenter).^2)/length(error(idxCenter)));
    rmseWhole = sqrt(sum(error(idxValid).^2)/length(error(idxValid)));
    disp(['rmse with gt mesh: Face:' num2str(rmseCenter) 'Head:' num2str(rmseWhole)]);
 
    %% Euclidian distance with ground truth mesh
%     euDist = mean(euclidian_dist(estShape(:,idxCenter),mesh(:,idxCenter)));
    %euDist = sqrt(sum((estShape(:,idxCenter) - mesh(:,idxCenter)).^2));
    error = sqrt(sum((estShape(:,idxCenter) - mesh(:,idxCenter)).^2));
    euDist = mean(error);
    
    figure(2);
    rp = defrp;
    rp.phi = 0; % frontal face
    rp.width = 500;
    rp.height = 500;
    rp.theta = 0.5*pi;
    rp.alpha = pi;
    shape = estShape(:);
%     error = euclidian_dist(gtShape, estShape,Segmentation.segValid);  
%     mae = mean(error);
    colormap jet;
    colorbar;
    cmap = colormap;
    maxerror = 10; % hottest err: 10 mm
    err = min(error/maxerror,1);
    coloridx=max(uint8(64*err),1);
    colors_N_3 = 200*repmat([1 1 1],[size(estShape,2),1]);
    colors_N_3(idxCenter,:) = 255*cmap(coloridx,:);
    colors_3_N = colors_N_3';
    colors_3N_1 = colors_3_N(:);
    tex = colors_3N_1;
    %tex = 200 * repmat([1;1;1],[1,size(estShape,2)]);
    display_face(shape, tex, faces, rp);  
    c = colorbar;
    c.Position = [0.9 0.02 0.04 0.96];
    c.FontSize = 12;
    c.TickLabels = maxerror*0:10;
    title(['3D face:estimated -vdist =' num2str(rmseParaCenter)]);
    
    %% show 3D gt shape
    figure(3);
    tex = 200 * repmat([1;1;1],[1,size(mesh,2)]);
    
    shape = mesh(:);
    display_face(shape, tex, faces, rp);
    title('3D face:mesh');
    %% show 3D est shape
    figure(4);
    tex = 200 * repmat([1;1;1],[1,size(estShape,2)]);
    shape = estShape(:);
    display_face(shape, tex, faces, rp);
    title('3D face:estimation');
    
    figure(5);
    gtShape = Data.worldR*gtShape + repmat(Data.worldT,[1 size(gtShape,2)]);
    tex = 200 * repmat([1;1;1],[1,size(gtShape,2)]);
    shape = gtShape(:);
    display_face(shape, tex, faces, rp);
    title('3D face:ground truth');
       
end
