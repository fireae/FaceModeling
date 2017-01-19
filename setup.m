function options = setup( )
%%  @paths    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cwd=cd;
cwd(cwd=='\')='/';

options.workDir    = cwd;
options.slash      = '/'; %% For linux

%% @library paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isdeployed
    
    %training
    addpath([cwd '/source/train/']);
    
    %testing
    addpath([cwd '/source/test/']);
    
    %common
    addpath([cwd '/common/io/']);
    addpath([cwd '/common/align/']);
    addpath([cwd '/common/string/']);
    addpath([cwd '/common/manifold/']);
    addpath([cwd '/common/err/']);
    addpath([cwd '/common/graphic/']);
    addpath([cwd '/common/vec/']);
    addpath([cwd '/common/desc/']);
    addpath([cwd '/common/regression/']);
    addpath([cwd '/common/flip/']);
    addpath([cwd '/common/eval/']);
    addpath([cwd '/common/BFM/']);
    %face detection algorithm from face++
    addpath([cwd '/facedetect/']);
    addpath([cwd '/model/']);
    
    %libs
    addpath([cwd '/lib/liblinear-2.1/matlab']);
    addpath([cwd '/lib/vlfeat/toolbox/']);
    addpath([cwd '/lib/vlfeat/bin/glnx86/']);
    
end

vl_setup();

%%  @data preparation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.datasetName = 'BFM';  %% 'lfpw', 'helen' or 'w300'

options.trainingImageDataPath = './data/BFM/trainfrontal/';
options.trainingTruthDataPath = './data/BFM/trainfrontal/';
                                   
options.testingImageDataPath  = './data/real/';
options.testingTruthDataPath  = './data/real/';

options.ResultDataPath = './Result/imgs/';
options.ResultFigurePath = './Result/cum/';
options.learningShape     = 0;
options.learningVariation = 0;
options.useBoundary = 0;


options.flipFlag          = 0;   % the flag of flipping

%%  @other folders  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.tempPath      = 'temp';
options.modelPath     = 'model';

%%  @cascade regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.n_cascades     = 10;         % the number of cascades

options.n_init_randoms = 2;         % the number of initial randoms

options.n_init_randoms_test = 1;    % the number of initial randoms
            
%%  @training configuration  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.canvasSize  = [400 400];
options.scaleFactor = 1;
%options.lambda     = [1 2 3 4 5 6 7] * 100;
%options.lambda      = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1] * 0.005;
options.lambda      = 0.005 * ones(1, options.n_cascades);
options.paraSize = 100;
%%  @feature configuration  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.descType  = 'hog'; % raw, hog (124 dims), dense_sift (128 dims)
options.descRawWin = 10;
options.descSize  = [20 20 20 20 20 20 20];
%options.descScale = [0.16 0.16 0.16 0.16 0.16 0.16 0.16 0.16 0.16 0.16];
options.descScale = 0.16 * ones(1,options.n_cascades);
options.descBins  =  4;

            

%%  @evaluation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%options.pts_eval = [1 3 2 4 9 17 11 12 18 10 19 20 21 23 24 25 28];
options.pts_eval = 1:68;
%options.pts_eval = 1:51;
%options.pts_eval = 1:29;

%% 300-W dataset (68)
% options.inter_ocular_left_eye_idx  = 37:42;
% options.inter_ocular_right_eye_idx = 43:48;

%% 300-W dataset (51)
%options.inter_ocular_left_eye_idx  = 20:25;
%options.inter_ocular_right_eye_idx = 26:31;

%options.inter_ocular_left_eye_idx  = 20;
%options.inter_ocular_right_eye_idx = 29;

%% LFPW dataset
%options.inter_ocular_left_eye_idx  = [9 11 13 14 17];
%options.inter_ocular_right_eye_idx = [10 12 15 16 18];

