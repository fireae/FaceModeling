function run_training ( )
close all;
tic;
%% loading the setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = setup( );

%[BFMmodel,~] = load_model(); %load BFM 3d face model
exist faces;
if ans == 0
    load('model/FWM5.mat');
end
%load  innerKeypointIndices.mat %load indices of inner key points on 3d face model
%boundIdxSet =[];


%% learn cascaded regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imgDir = options.trainingImageDataPath;
%ptsDir = options.trainingTruthDataPath;

%% loading data
disp('Loading training data...');
trainData = load_all_bd_data([imgDir 'bs000/'], options);
for id = 1:99
    path = [imgDir 'bs' sprintf('%03d',id) '/'];
    tmp = load_all_bd_data(path, options);
    trainData = [trainData;tmp];
end
clear tmp;
n_cascades = options.n_cascades;
LearnedCascadedModel{n_cascades}.R = [];
rms = zeros(n_cascades,1);

bv = [];
if options.useBoundary ==1
    %boundIdxSet = Landmarks.boundary;
    bv = zeros(length(trainData)*options.n_init_randoms, size(Landmarks.boundary,1));
end

for icascade = 1 : n_cascades
    
    options.current_cascade = icascade;
    
    %% learning single regressors
    if icascade == 1
        
        new_init_para = [];
        [R, new_init_para, rms(icascade), bv] = learn_single_regressor( ...
            Bilinear, faces, Landmarks, bv, trainData, new_init_para, options );
        LearnedCascadedModel{icascade}.R = R;   
        %% save other parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        LearnedCascadedModel{icascade}.n_cascades = n_cascades;
        LearnedCascadedModel{icascade}.descSize   = options.descSize;
        LearnedCascadedModel{icascade}.descBins   = options.descBins;
             
    else
        
        [R,new_init_para,rms(icascade), bv] = learn_single_regressor( ...
            Bilinear, faces, Landmarks, bv, trainData, new_init_para, options );     
        LearnedCascadedModel{icascade}.R = R;
        
    end   
    save('Result/Trained_RMS.mat' , 'rms');
end

%save('Result/Trained_RMS.mat' , 'rms');

save([options.modelPath options.slash ...
    'LearnedCascadedModel.mat'],'LearnedCascadedModel','options');
toc;
clear;
