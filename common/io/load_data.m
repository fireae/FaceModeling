function [data] = load_data( imgDir, dataDir , options )

slash = options.slash;

plist = dir([imgDir slash 'im*.*g']);
nlist = length(plist);
%nlist = 10;

data(nlist).shape   = [];
data(nlist).img     = '';

for ilist = 1 : nlist
    
    img_name = plist(ilist).name;
    %img_path = [imgDir slash img_name];
    data(ilist).img    = [];
    %load shape parameters
    data_name = [img_name(3:end-4) 'Para.data'];
    data_path = [dataDir slash data_name]; 
    fid = fopen(data_path,'rb');
    data(ilist).para = fread(fid,inf,'float');
    fclose(fid);
    %load transformation parameters
    data_name = [img_name(3:end - 4) 'Trans.data'];
    data_path = [dataDir slash data_name];
    fid = fopen(data_path,'rb');
    data(ilist).trans = fread(fid,inf,'float');
    fclose(fid);
      
    if 0
        img = imread(img_path);
         figure(1); imshow(img); hold on;
        hold off;
        pause;
    end
        
end

end
