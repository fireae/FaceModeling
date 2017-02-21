function shape = face_alignment( BFMmodel,keypoints,BoundIdxSet,...
    LearnedCascadedModel, Data, img, gt_trans, options,idata,cx,cy )



%% setup the fixed parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nData = length(Data);
n_init_randoms_test = options.n_init_randoms_test;
paraDim  = options.paraSize;

aligned_shape = zeros(n_init_randoms_test,paraDim);

%% randomize which shape is used for initial position
rIdx = randi([1,nData],n_init_randoms_test);

%% iterations of n initial points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ir = 1 : n_init_randoms_test
    
    %% get random positions and inital shape indexs
    idx    = rIdx(ir);
    init_para = Data(idx).para; %% get randomly shape from others
    init_trans = gt_trans;
    
    
    %% detect landmarks using cascaded regression
    aligned_shape(ir,:) = cascaded_regress( BFMmodel,keypoints,BoundIdxSet, ...
        LearnedCascadedModel, img, init_para, init_trans, options,idata,cx,cy );
    
    
    
end

if n_init_randoms_test == 1
    shape = aligned_shape;
else
    shape =aligned_shape;
end

end
