function bIndx = load_bIndex()

for bv = 1:17
    fid = fopen(['boundary' num2str(bv-1) '.data'],'rb');
    buffer = fread(fid,inf,'int');
    fclose(fid);
    buffer = buffer(buffer ~=0);
    buffer = buffer+ones(size(buffer));
    
    bIndx{bv} = buffer;
end
bIndx{1} = flipud(bIndx{1});
bIndx{2} = flipud(bIndx{2});
bIndx{7} = flipud(bIndx{7});
bIndx{8} = flipud(bIndx{8});
end