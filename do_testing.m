function do_testing ( )
%% clear all;
clc;close all;
%% loading the setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = setup( );
%% loading FWM model
exist faces;
if ans == 0
    load('model/FWM5.mat');
end
%% loading Learned Model
load( [options.modelPath options.slash 'LearnedCascadedModel.mat'] );
%% loading training paras for randomly initialize shapes.
imgTrainDir = options.trainingImageDataPath;
% TrainingData = load_all_bd_data([imgTrainDir 'bs000/'], options);
% for id = 1:99
%     path = [imgTrainDir 'bs' sprintf('%03d',id) '/'];
%     tmp = load_all_bd_data(path, options);
%     TrainingData = [TrainingData;tmp];
% end
% clear tmp;
%% load testing data
imgDir = options.testingImageDataPath;
Data = load_all_bd_data([imgDir 'bs100/'],options);
for id = 101:104
    path = [imgTrainDir 'bs' sprintf('%03d',id) '/'];
    tmp = load_all_bd_data(path, options);
    Data = [Data;tmp];
end
clear tmp;
nData = length(Data);

%% evaluating on whole data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pdc = zeros(nData,1);
rmse = zeros(nData,1);
estimation = cell(1,nData);
CrTensor = tensor(Bilinear.Cr);

for idata = 1 : nData
    
    disp(['Image: ' num2str(idata)]);
    
    %% estimate parameters based on the learning regressor
    [estID,R,T,s] = face_alignment( Bilinear, faces, Landmarks, ...
        LearnedCascadedModel, Data{idata}, idata, options);
    estimation{idata}.ID = estID;
    estimation{idata}.R = R;
    estimation{idata}.T = T;
    estimation{idata}.s = s;
    save('Result/estimation.mat','estimation');    
   [error1,error2] = visualize_result(Data{idata},estimation{idata},CrTensor, faces, Landmarks, Segmentation, options);
   pdc(idata) = error1;
   rmse(idata) = error2;
    
end

%% displaying CED
if 1 
    %% save measure results
    save('Result/pdc.mat', 'pdc');
    save('Result/rmse.mat', 'rmse');

    %% plot projection distance
    figure(4);
    visualize_ced(pdc,options);
    figure(5);
    visualize_ced(rmse,options);
    
    %% displaying rms errors
    disp(['pdc average: ' num2str(mean(pdc))]);
    disp(['rmse average: ' num2str(mean(rmse))]);
end


