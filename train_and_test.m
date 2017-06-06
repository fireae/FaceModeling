function train_and_test ( )
%%%%%%%%%%%%%%%%%%%%%%%%% train
close all;
tic;
%% loading the setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = setup();
parpool('local',40);
% load FaceWarehouse Model
if (exist('faces','var')==0) 
    tmp = load('model/FWM5.mat');
    faces = tmp.faces;
    Bilinear = tmp.Bilinear;
    Landmarks = tmp.Landmarks;
    Segmentation = tmp.Segmentation;
    clear tmp;
end

%% learn cascaded regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imgDir = options.trainingImageDataPath;
%ptsDir = options.trainingTruthDataPath;

%% loading data 
disp('Loading training data...');
trainData = load_bd_data([imgDir 'bs000/'], options);
parfor id = 1:29
    path = [imgDir 'bs' sprintf('%03d',id) '/'];
    tmp = load_bd_data(path, options);
    trainData = [trainData;tmp];
end
parfor id = 35:104
    path = [imgDir 'bs' sprintf('%03d',id) '/'];
    tmp = load_bd_data(path, options);
    trainData = [trainData;tmp];
end
clear tmp;
n_cascades = options.n_cascades;
LearnedCascadedModel{n_cascades}.Regression = [];
rms = zeros(n_cascades,1);


for icascade = 1 : n_cascades
    
    options.current_cascade = icascade;
    
    %% learning single regressors
    if icascade == 1
        
        newInitPara = [];
        newInitTran = [];
        [Regression, newInitPara, newInitTran, rms(icascade)] = learn_single_regressor( ...
            Bilinear, faces, Landmarks, trainData, newInitPara, newInitTran, options );
        LearnedCascadedModel{icascade}.Regression = Regression;   
        %% save other parameters 
        LearnedCascadedModel{icascade}.n_cascades = n_cascades;
        LearnedCascadedModel{icascade}.descSize   = options.descSize;
        LearnedCascadedModel{icascade}.descBins   = options.descBins;
             
    else
        
        [Regression,newInitPara,newInitTran, rms(icascade)] = learn_single_regressor( ...
            Bilinear, faces, Landmarks, trainData, newInitPara, newInitTran, options );     
        LearnedCascadedModel{icascade}.Regression = Regression;
        
    end   
    
end

save([options.ResultPath 'TrainLoss.mat'] , 'rms');

save([options.ResultPath 'LearnedCascadedModel.mat'],'LearnedCascadedModel');
toc;

% %%%%%%%%%%%%%%%%%%%%%%%%%validate training data
% load([options.ResultPath 'LearnedCascadedModel.mat']);
% nData = length(trainData);
% pdc = zeros(nData,1);
% rmse1center = zeros(nData,1);
% rmse1whole = zeros(nData,1);

% validation = cell(1,nData);
% CrTensor = tensor(Bilinear.Cr);

% %load([options.ResultPath 'validation.mat']);

% for idata = 1:nData%nData
    
%     disp(['Image: ' num2str(idata)]);
    
%     %% estimate parameters based on the learning regressor
%     [estID, estEP, R, T, s] = face_alignment( Bilinear, faces, Landmarks, ...
%        LearnedCascadedModel, trainData{idata}, idata, options);    
%     validation{idata}.EP = estEP;
%     validation{idata}.ID = estID;
%     validation{idata}.R = R;
%     validation{idata}.T = T;
%     validation{idata}.s = s;
    
%     if trainData{idata}.flag ==1
%     [error1,error2,error3] = validate_data(trainData{idata},validation{idata},CrTensor, faces, Landmarks, Segmentation, options);
%     pdc(idata) = error1;
%     rmse1center(idata) = error2;
%     rmse1whole(idata) = error3;
%     end
% %     pause;
    
% end

% save([options.ResultPath 'validation.mat'],'validation');
% disp(['pdc average: ' num2str(mean(pdc(pdc>0)))]);
% disp(['rmse1center average: ' num2str(mean(rmse1center(rmse1center>0)))]);
% disp(['rmse1whole average: ' num2str(mean(rmse1whole(rmse1whole>0)))]);

% %% save measure results
% save([options.ResultPath 'valid_pdc.mat'], 'pdc');
% save([options.ResultPath 'valid_rmse1center.mat'], 'rmse1center');
% save([options.ResultPath 'valid_rmse1whole.mat'], 'rmse1whole');

% clear;
%%%%%%%%%%%%%%%%%%%%%%%%% test
clc;close all;
% loading the setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = setup( );
% loading FWM model
if (exist('faces','var')==0) 
    tmp = load('model/FWM5.mat');
    faces = tmp.faces;
    Bilinear = tmp.Bilinear;
    Landmarks = tmp.Landmarks;
    Segmentation = tmp.Segmentation;
    clear tmp;
end
%% loading Learned Model
load([options.ResultPath 'LearnedCascadedModel.mat']);
LearnedCascadedModel = LearnedCascadedModel;
%% load testing data
imgDir = options.testingImageDataPath;
Data = load_all_bd_data([imgDir 'bs030/'],options);
parfor id = 30:34
    path = [imgDir 'bs' sprintf('%03d',id) '/'];
    tmp = load_all_bd_data(path, options);
    Data = [Data;tmp];
end
clear tmp;
% Data = load_data('./data/real/', options);
nData = length(Data);

%% evaluating on whole data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pdc = zeros(nData,1);
rmse1center = zeros(nData,1);
rmse1whole = zeros(nData,1);
rmse2center = zeros(nData,1);
rmse2whole = zeros(nData,1);
eudist = zeros(nData,1);
estimation = cell(1,nData);
CrTensor = tensor(Bilinear.Cr);
%load([options.ResultPath 'estimation.mat']);
parfor idata = 1:nData%nData
%     idata = randi([1,nData],1);
    disp(['Image: ' num2str(idata)]);
    
    % estimate parameters based on the learning regressor
    [estID, estEP, R, T, s] = face_alignment( Bilinear, faces, Landmarks, ...
       LearnedCascadedModel, Data{idata}, idata, options);    
    estimation{idata}.EP = estEP;
    estimation{idata}.ID = estID;
    estimation{idata}.R = R;
    estimation{idata}.T = T;
    estimation{idata}.s = s;
    
%     img = Data{idata}.imgCrop;
% %     allKey = Landmarks.inner; 
%     estShape = coeff2shape(CrTensor, estID, estEP); 
%     transformedestShape = R*estShape+repmat(T,[1,size(estShape,2)]);
%     tmp= get_boundary_vertex(transformedestShape, faces, Landmarks.boundary);
%     allKey = [Landmarks.inner tmp];
%     keyPoints = estShape(:,allKey); 
%     estProj = cal_weak_perspective(keyPoints, s, R, T);
%     figure(2);
%     imshow(img);
%     hold on;plot(estProj(1, :),estProj(2,:),'r*');
% %     hold on; plot(gtProj(1,:),gtProj(2,:),'g*');
%     hold on;plot(Data{idata}.lm2d(:,1),Data{idata}.lm2d(:,2),'b*');
%     hold off;
%     saveas(gcf,[options.ResultPath num2str(idata) '_proj.jpg']);
%     
%     rp = defrp;
%     rp.phi = 0; % frontal face
%     rp.width = 500;
%     rp.height = 500;
%     rp.theta = 0.5*pi;
%     rp.alpha = pi;
%     figure(3);
%     tex = 200 * repmat([1;1;1],[1,size(estShape,2)]);
%     shape = estShape(:);
%     display_face(shape, tex, faces, rp);
%     title('3D face:estimation');
%     saveas(gcf,[options.ResultPath num2str(idata) '_shape.jpg']);

%     if Data{idata}.flag ==1
    [error1,error2,error3,error4,error5,error6] = visualize_result(Data{idata},estimation{idata},CrTensor, faces, Landmarks, Segmentation, options);
     saveas(figure(1),[options.ResultPath 'imgs\' num2str(idata) '_img.jpg']);
     saveas(figure(2),[options.ResultPath 'imgs\' num2str(idata) '_errormap.jpg']);
     saveas(figure(3),[options.ResultPath 'imgs\' num2str(idata) '_mesh.jpg']);
     saveas(figure(4),[options.ResultPath 'imgs\' num2str(idata) '_estshape.jpg']);
     saveas(figure(5),[options.ResultPath 'imgs\' num2str(idata) '_gtshape.jpg']);
    pdc(idata) = error1;
    rmse1center(idata) = error2;
    rmse1whole(idata) = error3;
    rmse2center(idata) = error4;
    rmse2whole(idata) = error5;
    eudist(idata) = error6;
%     end
%     pause;
    
end
save([options.ResultPath 'estimation.mat'],'estimation');

    disp(['pdc average: ' num2str(mean(pdc(pdc>0)))]);
    disp(['rmse1center average: ' num2str(mean(rmse1center(rmse1center>0)))]);
    disp(['rmse1whole average: ' num2str(mean(rmse1whole(rmse1whole>0)))]);
    disp(['rmse2center average: ' num2str(mean(rmse2center(rmse2center>0)))]);
    disp(['rmse2whole average: ' num2str(mean(rmse2whole(rmse2whole>0)))]);
    disp(['Euclidian Distance:' num2str(mean(eudist(eudist>0)))]);

    %% save measure results
    save([options.ResultPath 'pdc.mat'], 'pdc');
    save([options.ResultPath 'rmse1center.mat'], 'rmse1center');
    save([options.ResultPath 'rmse1whole.mat'], 'rmse1whole');
    save([options.ResultPath 'rmse2center.mat'], 'rmse2center');
    save([options.ResultPath 'rmse2whole.mat'], 'rmse2whole');
    save([options.ResultPath 'eudist.mat'], 'eudist');
% %% displaying CED    
% if 0 
%     %% plot projection distance
%     figure(5);
%     visualize_ced(pdc,options);
%     saveas(figure(4),[options.ResultPath 'pdc.jpg']);
%     figure(6);
%     visualize_ced(rmse,options);
%     saveas(figure(5),[options.ResultPath 'rmse.jpg']);
    
%     %% displaying rms errors
%     disp(['pdc average: ' num2str(mean(pdc))]);
%     disp(['rmse average: ' num2str(mean(rmse))]);
% end
