function display_face (shp, tex, tl, rp)
	
	shp = reshape(shp, [ 3 prod(size(shp))/3 ])'; 
	tex = reshape(tex, [ 3 prod(size(tex))/3 ])'; 
    %% Only show facial part
     load nonface3.mat
     tex(idx,:) = nan;
     shp(idx,:) = nan;
    %% only show facial part
	tex = min(tex, 255);
    
    
%     BoundIdxSet = load_bIndex();
%     Bound = zeros(17,1);
%     tmpidx = 20 * ones(17,1);
%     tmpidx(3) = 25;
%     tmpidx(4) = 35;
%     tmpidx(5) = 45;
%     tmpidx(6) = 35;
%     
%     tmpidx(8:10) = 10;
%     for i = 1:size(BoundIdxSet,2)
%         tmp = BoundIdxSet{i};
%         %tex(tmp(1),:) = repmat([255 0 0],size(BoundIdxSet{i},1),1);
%         tex(tmp(tmpidx(i)),:) = [255 0 0];
%         Bound(i) = tmp(tmpidx(i));
%     end
%     
%     
%     tmp = BoundIdxSet{5};
   % tex(tmp(23:end),:) = repmat([255 0 0],size(BoundIdxSet{5},1)-22,1);
   


	set(gcf, 'Renderer', 'opengl');
	fig_pos = get(gcf, 'Position');
	fig_pos(3) = rp.width;
	fig_pos(4) = rp.height;

	mesh_h = trimesh(...
		tl, shp(:, 1), shp(:, 3), shp(:, 2), ...
		'EdgeColor', 'none', ...
		'FaceVertexCData', tex/255, 'FaceColor', 'interp', ...
		'FaceLighting', 'phong' ...
	);

	set(gca, ...
        'PlotBoxAspectRatio', [ 1 1 1 ], ...
		'DataAspectRatio', [1 1 1 ], ...
		'Units', 'pixels', ...
		'GridLineStyle', 'none', ...
        'Position', [ 0 0 fig_pos(3) fig_pos(4) ], ...
		'Visible', 'on', 'box', 'on', ...
		'Projection', 'perspective', ...
        'CameraPosition',[0 0 0]...
        ); 
%         ...
% 		
		
	
	set(gcf, 'Color', [ 0 0 0 ]);
    out = get(gca,'CameraPosition');
    out = get(gca,'CameraTarget');
    set(gca,'CameraTarget',[0 0 0]);
	view(-180, 0);

%     set(gca,'CameraViewAngleMode','manual');

	material([.5, .5, .1 1  ])
	camlight('headlight');

	
%% ------------------------------------------------------------CALLBACK--------
function resizeCallback (obj, eventdata)
	
	fig = gcbf;
	fig_pos = get(fig, 'Position');

	axis = findobj(get(fig, 'Children'), 'Tag', 'Axis.Head');
	set(axis, 'Position', [ 0 0 fig_pos(3) fig_pos(4) ]);
	
