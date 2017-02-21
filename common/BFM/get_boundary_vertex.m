function bIdx = get_boundary_vertex(shape, faces, BoundIdxSet)
%bIdx = zeros(1,size(BoundIdxSet,2));
%otherpoints = [];
for i = 1:size(BoundIdxSet,2)
    CurrentBset = BoundIdxSet{i};
    norm = zeros(3,size(CurrentBset,1));
    for j = 1:size(CurrentBset,1)
        vIdx = CurrentBset(j);
        norm(:,j) = cal_vertex_norm(shape, faces, vIdx);
    end
    [~,currentIdx] = min(abs(norm(3,:)));
    bIdx(i) = CurrentBset(currentIdx);
    %bIdx(i) = CurrentBset(1);
%     CurrentBset(currentIdx) = [];
%     otherpoints = cat(2,otherpoints,CurrentBset');
end

%%use fixed boundary indices
% load FixBoundIdx
% bIdx = Bound';

end