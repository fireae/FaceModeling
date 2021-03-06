function bIdx = get_boundary_vertex(shape, faces, BoundIdxSet)
%bIdx = zeros(1,size(BoundIdxSet,2));
%otherpoints = [];
%BoundaIdxSet [N K] N:number of landmarks K:number of possible keypoints
%for one landmarks
bIdx = zeros(1,size(BoundIdxSet,1));
for i = 1:size(BoundIdxSet,1)
    Currentset = BoundIdxSet(i,:);
    norm = zeros(3,size(Currentset,2));
    for j = 1:size(Currentset,2)
        vIdx = Currentset(j);
        norm(:,j) = cal_vertex_norm(shape, faces, vIdx);
    end
    [~,currentIdx] = min(abs(norm(3,:)));
    bIdx(i) = Currentset(currentIdx);
    %bIdx(i) = CurrentBset(1);
%     CurrentBset(currentIdx) = [];
%     otherpoints = cat(2,otherpoints,CurrentBset');
end

%%use fixed boundary indices
% load FixBoundIdx
% bIdx = Bound';

end