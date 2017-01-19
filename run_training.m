function run_training ( )
tic;
%% loading the setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = setup( );

[BFMmodel,~] = load_model(); %load BFM 3d face model
load  innerKeypointIndices.mat %load indices of inner key points on 3d face model
boundIdxSet =[];


%% learn cascaded regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imgDir = options.trainingImageDataPath;
ptsDir = options.trainingTruthDataPath;

%% loading data
disp('Loading training data...');
trainData = load_all_data2(imgDir, ptsDir, options );
n_cascades = options.n_cascades;
LearnedCascadedModel{n_cascades}.R = [];
rms = zeros(n_cascades,1);

bv = [];
if options.useBoundary ==1
    boundIdxSet = load_bIndex();
    bv = zeros(length(trainData)*options.n_init_randoms, size(boundIdxSet,2));
end

for icascade = 1 : n_cascades
    
    options.current_cascade = icascade;
    
    %% learning single regressors
    if icascade == 1
        
        new_init_para = [];
        [R, new_init_para, rms(icascade), bv] = learn_single_regressor( ...
            BFMmodel,  innerKeypointIndices, boundIdxSet, bv, trainData, new_init_para, options );
        LearnedCascadedModel{icascade}.R = R;   
        %% save other parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        LearnedCascadedModel{icascade}.n_cascades = n_cascades;
        LearnedCascadedModel{icascade}.descSize   = options.descSize;
        LearnedCascadedModel{icascade}.descBins   = options.descBins;
             
    else
        
        [R,new_init_para,rms(icascade), bv] = learn_single_regressor( ...
            BFMmodel,  innerKeypointIndices, boundIdxSet, bv, trainData, new_init_para, options );     
        LearnedCascadedModel{icascade}.R = R;
        
    end   
    save('Result/Trained_RMS.mat' , 'rms');
end

%save('Result/Trained_RMS.mat' , 'rms');

save([options.modelPath options.slash ...
    'LearnedCascadedModel.mat'],'LearnedCascadedModel');
toc;
clear;
