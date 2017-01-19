#function do_testing ( )

%clear all;

%% loading the setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = setup( );
%% loading training data

[BFMmodel,msz] = load_model();
load  innerKeypointIndices.mat
BoundIdxSet =[];
if options.useBoundary ==1
    BoundIdxSet = load_bIndex();
end
load( [options.modelPath options.slash 'LearnedCascadedModel.mat'] );

%% loading training paras for randomly initialize shapes.
imgTrainDir = options.trainingImageDataPath;
ptsTrainDir = options.trainingTruthDataPath;
TrainingData = load_data( imgTrainDir, ptsTrainDir, options );
%% load testing data
imgDir = options.testingImageDataPath;
ptsDir = options.testingTruthDataPath;
Data  = load_real_data(imgDir, ptsDir, options );
nData = length(Data);
%nData = 20;
%% evaluating on whole data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%err = zeros(nData,1);
dist = zeros(nData,1);
verror = zeros(nData,1);
for idata = 1 : nData
    close all;
    disp(['Image: ' num2str(idata)]);
    
    %% information of one image
    img        = Data{idata}.img_gray;
    
    truePara = Data{idata}.para_gt;%gt:[207,1];
    trueTrans = Data{idata}.trans_gt;
    imsize = size(img);
    
    %% estimate parameters based on the learning regressor
    estPara = face_alignment( BFMmodel, innerKeypointIndices, BoundIdxSet, ...
        LearnedCascadedModel, TrainingData, img, trueTrans, options, idata,msz );
    
    allKey = innerKeypointIndices;    
    
    
    allKey = innerKeypointIndices;          
    estShape = Reconstruct_face(BFMmodel, estPara', trueTrans);
    if options.useBoundary ==1
        bv = get_boundary_vertex(estShape, BFMmodel.tl, BoundIdxSet);
        allKey = [innerKeypointIndices bv];
    end    
    estProj = CalKeyProj(estShape, trueTrans, allKey, imsize); 
    
    if (~isnan(truePara))
    gtShape = Reconstruct_face(BFMmodel, truePara, trueTrans);
    if options.useBoundary ==1
        bv = get_boundary_vertex(gtShape, BFMmodel.tl, BoundIdxSet);
        allKey = [innerKeypointIndices bv];
    end

    trueProj = CalKeyProj(gtShape, trueTrans, allKey,imsize);
    dist(idata) = norm_Edist(trueProj, estProj, options);
    verror(idata) = vdc(gtShape, estShape);  
    end
    figure(1); imshow(img); hold on;
        
        plot(estProj(:,1),estProj(:,2),'g.');
        
        beta  = zeros(msz.n_tex_dim, 1); % use mean texture for rendering
        tex    = coef2object( beta,  BFMmodel.texMU,   BFMmodel.texPC,  BFMmodel.texEV );       
        tmpshape = coef2object(estPara', BFMmodel.shapeMU, BFMmodel.shapePC, BFMmodel.shapeEV);
        
        rp     = defrp;
        rp.phi = 0.5;
        rp.dir_light.dir = [0;1;1];
        rp.dir_light.intens = 0.6*ones(3,1);
        figure(2); 
        display_face(tmpshape, tex, BFMmodel.tl, rp);
        title('3D face:estimated');
        saveas(gcf,[options.ResultDataPath num2str(idata) 'estimate3D.jpg']);
    if 0
        
        figure(1); imshow(img); hold on;
        plot(trueProj(:,1),trueProj(:,2),'r.');hold on;
        
        plot(estProj(:,1),estProj(:,2),'g.');
        hold off;
        title(['Key points projection Distance:' num2str(dist(idata))]);
        saveas(gcf,[options.ResultDataPath num2str(idata) 'proj.jpg']);
        legend('ground truth','estimate');
        
        % Render it
        beta  = zeros(msz.n_tex_dim, 1); % use mean texture for rendering
        tex    = coef2object( beta,  BFMmodel.texMU,   BFMmodel.texPC,  BFMmodel.texEV );       
        tmpshape = coef2object(estPara', BFMmodel.shapeMU, BFMmodel.shapePC, BFMmodel.shapeEV);
        
        rp     = defrp;
        rp.phi = 0.5;
        rp.dir_light.dir = [0;1;1];
        rp.dir_light.intens = 0.6*ones(3,1);
        figure(2); 
        display_face(tmpshape, tex, BFMmodel.tl, rp);
        title('3D face:estimated');
        saveas(gcf,[options.ResultDataPath num2str(idata) 'estimate3D.jpg']);
        
        tmpshape = coef2object(truePara, BFMmodel.shapeMU, BFMmodel.shapePC, BFMmodel.shapeEV);
        figure(3); 
        display_face(tmpshape, tex, BFMmodel.tl, rp);
        title('3D face:ground truth');
        saveas(gcf,[options.ResultDataPath num2str(idata) 'ground3D.jpg'])
        %pause;
    end  
end

save('Result/dist.mat', 'dist');
save('Result/verror.mat', 'verror');

%% displaying CED


if 1 
    x = 0:0.05:1;
x = x*(max(dist) - min(dist))+min(dist);
CED = zeros(length(x),1);
c = 0;

    for thres = x
        c = c + 1;
        idx = find(dist <= thres);
        CED(c) = length(idx)/nData;
    end

    figure(5);
    plot( x, CED, 'LineWidth', 2 , 'MarkerEdgeColor','r');
    axis([min(dist) max(dist) 0 1])
    title(['CED: Dist average: ' num2str(mean(dist))]);
    grid on;
    saveas(gcf,[options.ResultFigurePath 'CED.jpg']);

    %%show vertex distance curve
    x = 0:0.05:1;
    x = x * (max(verror) - min(verror)) + min(verror);
    V = zeros(length(x),1);
    c = 0;

    for thres = x
        c = c + 1;
        idx = find(verror <= thres);
        V(c) = length(idx)/nData;
    end

    figure(6);
    plot( x, V, 'LineWidth', 2 , 'MarkerEdgeColor','r');
    axis([min(verror) max(verror) 0 1])
    title(['VDC: Vetex distance error: ' num2str(mean(verror))]);
    grid on;
    saveas(gcf,[options.ResultFigurePath 'VDC.jpg']);
end
%% displaying rms errors
disp(['Dist average: ' num2str(mean(dist))]);
disp(['VDist average: ' num2str(mean(verror))]);

