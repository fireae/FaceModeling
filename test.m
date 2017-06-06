%close all;
%git test
oldlist = dir(['.\data\BFM\trainfrontal\' '*.jpg']);
for i = 1:length(oldlist)
    im = imread(['.\data\BFM\traincrop\im' num2str(i) '.jpg']);
    imwrite(im,['.\data\BFM\tmp\' oldlist(i).name])
    
    load (['.\data\BFM\traincrop\cam' num2str(i) '.mat'])
    camera.cx = camera.cx - 320;
    camera.cy = camera.cy - 240;
    save(['.\data\BFM\tmp\cam' oldlist(i).name(3:end - 4) '.mat'],'camera');
end


% options = setup( );
% %[BFMmodel,~] = load_model(); %load BFM 3d face model
% %% learn cascaded regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% imgDir = './data/real/';
% ptsDir = './data/real/';
% 
% %% loading data
% disp('Loading training data...');
% Data = load_real_data(imgDir, ptsDir, options);
% for i = 1:length(Data)
% %     [bbox,~] = facedetect(Data{i}.img);
% %     cropim = imcrop(Data{i}.img,bbox);
% %     imwrite(cropim, ['data\real\crop\' num2str(i) 'im.jpg']);
% load (['data\real\crop\' num2str(i) 'cam.mat']);
% 
%     camera.cx = size(Data{i}.img,2)/2.0 - (Data{i}.centery - camera.cx);  %cx: center along horizontal direction
%     camera.cy = size(Data{i}.img,1)/2.0 - (Data{i}.centerx - camera.cy);  %cy: center along vertical direction
%     %trans = Data{i}.trans_gt;
%     %para = Data{i}.para_gt;
%     save (['data\real\crop\' num2str(i) 'cam.mat'],'camera');
% %     save (['data\real\crop\' num2str(i) 'Trans.mat'],'trans');
% %     save (['data\real\crop\' num2str(i) 'Para.mat'],'para');
%     
%     
%     clear camera;
%     clear trans;
%     clear para;
% end
% load  innerKeypointIndices.mat
% 
% image = imread('bs000_LFAU_12LW_0.png');
% load bs000_LFAU_12LW_0_campar.mat
% 
% para = zeros(100,1);
% 
% trans = zeros(7,1);
% 
% trans(1:3) = Rot2Ang(Rc).*(180/pi); %transfer rotation matrix to angle degree
% tmp = Ang2Rot(trans);
% trans(4:6) = Tc;
% trans(7) = fx;
% 
% shape = Reconstruct_face(BFMmodel,para,trans);
% 
% proj = CalKeyProj(shape,trans,innerKeypointIndices,size(image),cx,cy);
% 
% figure; imshow(image); 
% hold on;
% plot(proj(:,1),proj(:,2),'r.');
% hold off;
% tmp = Ang2Rot(trans.*(180/pi));



% load Result12/test/dist.mat
% load Result12/test/verror.mat
% dist1 = dist;
% verror1 = verror;
% clear dist
% clear verror
% 
% load Result12/train/dist.mat
% load Result12/train/verror.mat
% dist1t = dist;
% verror1t = verror;
% clear dist
% clear verror
% 
% load Result15/test/dist.mat
% load Result15/test/verror.mat
% dist2 = dist;
% verror2 = verror;
% clear dist
% clear verror
% 
% load Result15/train/dist.mat
% load Result15/train/verror.mat
% dist2t = dist;
% verror2t = verror;
% clear dist
% clear verror
% 
% % load Result18/test/dist.mat
% % load Result18/test/verror.mat
% % dist3 = dist;
% % verror3 = verror;
% % clear dist
% % clear verror
% 
% % dist = dist(1:200);
% % verror = verror(1:200);
% nData = length(dist1);
% x = 0:0.06:1;
% x = x*(max([dist1;dist2]) - min([dist1;dist2]))+min([dist1;dist2]);
% CED1 = zeros(length(x),1);
% CED2 = zeros(length(x),1);
% CED1t = zeros(length(x),1);
% CED2t = zeros(length(x),1);
% CED3 = zeros(length(x),1);
% c = 0;
% 
% for thres = x
%         c = c + 1;
%         idx = find(dist1 <= thres);
%         CED1(c) = length(idx)/nData;
%         idx = find(dist1t <= thres);
%         CED1t(c) = length(idx)/nData;
%         idx = find(dist2 <= thres);
%         CED2(c) = length(idx)/nData;
%         idx = find(dist2t <= thres);
%         CED2t(c) = length(idx)/nData;
% %         idx = find(dist3 <= thres);
% %         CED3(c) = length(idx)/nData;
%     end
% 
%     figure(1);
%     plot( x, CED1, 'LineWidth', 2 , 'color','r');
%     axis([min([dist1;dist2]) max([dist1;dist2]) 0 1.1])
%     hold on;
%     plot( x, CED1t, 'LineStyle','- . ', 'LineWidth', 2,'color','r');
%     plot( x, CED2, 'LineWidth', 2 , 'color','b');
%     plot( x, CED2t, 'LineStyle','-.','LineWidth' ,2 ,'color','b');
% %     plot( x, CED3,'LineWidth' ,2 ,'color','k');
%     hold off;
%     title(['CED: Projection Dist Comparison' ]);
%     grid on;
%     legend('5Iter-test','5Iter-train','15Iter - test','15Iter - train','Location','SouthEast');
%     saveas(gcf,'CED.jpg');
% 
%     %%show vertex distance curve
% x = 0:0.06:1;
% x = x*(max([verror1;verror1t;verror2;verror2t]) - min([verror1;verror1t;verror2;verror2t]))+min([verror1;verror1t;verror2;verror2t]);
% CED1 = zeros(length(x),1);
% CED2 = zeros(length(x),1);
% CED1t = zeros(length(x),1);
% CED2t = zeros(length(x),1);
% c = 0;
% 
% for thres = x
%         c = c + 1;
%         idx = find(verror1 <= thres);
%         CED1(c) = length(idx)/nData;
%         idx = find(verror1t <= thres);
%         CED1t(c) = length(idx)/nData;
%         idx = find(verror2 <= thres);
%         CED2(c) = length(idx)/nData;
%         idx = find(verror2t <= thres);
%         CED2t(c) = length(idx)/nData;
% %         idx = find(verror3 <= thres);
% %         CED3(c) = length(idx)/nData;
%     end
% 
%     figure(2);
%     plot( x, CED1, 'LineWidth', 2 , 'color','r');
%     axis([min([verror1;verror1t;verror2;verror2t]) max([verror1;verror1t;verror2;verror2t]) 0 1.1])
%     hold on;
%     plot( x, CED1t, 'LineStyle','- . ', 'LineWidth', 2,'color','r');
%     plot( x, CED2, 'LineWidth', 2 , 'color','b');
%     plot( x, CED2t, 'LineStyle','-.','LineWidth' ,2 ,'color','b');
% %     plot( x, CED3, 'LineWidth', 2 , 'color','k');
%     hold off;
%     title(['CED:  Vertex Dist Comparison' ]);
%     grid on;
%     legend('5Iter-test','5Iter-train','15Iter - test','15Iter - train','Location','SouthEast');
%     saveas(gcf,'VED.jpg');
% % [BFMmodel,msz] = load_model();
% % % para = randn([207,1]);
% % % fid = fopen('test.data','rb');
% % % para = fread(fid,inf,'float');
% % % fclose(fid);
% % % 
% % % Recon_shape = Reconstruct_face(BFMmodel,para);
% % % BoundIdxSet = load_bIndex();
% % % bv = get_boundary_vertex(Recon_shape, BFMmodel.tl, BoundIdxSet);
% % 
% % alpha  = zeros(199, 1);
% % beta  = zeros(msz.n_tex_dim, 1);
% % shape  = coef2object(alpha, BFMmodel.shapeMU, BFMmodel.shapePC, BFMmodel.shapeEV );
% % tex    = coef2object( beta,  BFMmodel.texMU,   BFMmodel.texPC,  BFMmodel.texEV );
% % 
% % % Render it
% %         rp     = defrp;
% %         rp.phi = 0.5;
% %         rp.dir_light.dir = [0;1;1];
% %         rp.dir_light.intens = 0.6*ones(3,1);
% %         figure(2); 
% %         display_face(shape, tex, BFMmodel.tl, rp);
