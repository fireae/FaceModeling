function [bbox, landmarks ] = facedetect( img )
    %FACEDETECT Summary of this function goes here
    %   Detailed explanation goes here
    % compile.m should work for Linux and Mac.
    % To Windows users:
    % If you are using a Windows machine, please use the basic convolution (fconv.cc).
    % This can be done by commenting out line 13 and uncommenting line 15 in
    % compile.m
    %compile;
    % load and visualize model
    % Pre-trained model with 146 parts. Works best for faces larger than 80*80
    
    load face_p146_small.mat
    

    % % Pre-trained model with 99 parts. Works best for faces larger than 150*150
    % load face_p99.mat

    % % Pre-trained model with 1050 parts. Give best performance on localization, but very slow
    % load multipie_independent.mat

    % 5 levels for each octave
    model.interval = 5;
    % set up the threshold
    model.thresh = min(-0.65, model.thresh);

    % define the mapping from view-specific mixture id to viewpoint
    %     if length(model.components)==13 
    %         posemap = 90:-15:-90;
    %     elseif length(model.components)==18
    %         posemap = [90:-15:15 0 0 0 0 0 0 -15:-15:-90];
    %     else
    %         error('Can not recognize this model');
    %     end

    clf; axis image; axis off; drawnow;
    
    tic;
    bs = detect(img, model, model.thresh);
    
    bs = nms_face(bs,0.3);
    bs = clipboxes(img, bs);
    bs = bs(1);
    dettime = toc;
    
    bbox = zeros(1,4);
    bbox(1) = min(bs.xy(:,1));
    bbox(2) = min(bs.xy(:,2));
    bbox(3) = max(bs.xy(:,3)) - min(bs.xy(:,1));
    bbox(4) = max(bs.xy(:,4)) - min(bs.xy(:,2));
    landmarks = zeros(2,size(bs.xy,1));
    for i = size(bs.xy,1):-1:1;
        x1 = bs.xy(i,1);
        y1 = bs.xy(i,2);
        x2 = bs.xy(i,3);
        y2 = bs.xy(i,4);
        %line([x1 x1 x2 x2 x1]', [y1 y2 y2 y1 y1]', 'color', 'b', 'linewidth', 1);
        landmarks(1,i) = (x1+x2)/2;
        landmarks(2,i) = (y1+y2)/2;
        %plot((x1+x2)/2,(y1+y2)/2,'r.','markersize',15);
    end
    
    % show highest scoring one
    figure(1);imshow(img);
    %figure,showboxes(img, bs(1),posemap),title('Highest scoring detection');
    hold on;
    rectangle('Position',bbox);
    hold on;
    plot(landmarks(1,:),landmarks(2,:),'r.');
    % show all
    %figure,showboxes(im, bs,posemap),title('All detections above the threshold');
    
    fprintf('Detection took %.1f seconds\n',dettime);
end

