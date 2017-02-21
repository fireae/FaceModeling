function [BFMmodel msz] = load_model()
 % global BFMmodel;
  %if isempty(BFMmodel);
    BFMmodel = load('01_MorphableModel.mat');
  %end
  BFMmodel.shapeMU = BFMmodel.shapeMU/1000;
  msz.n_shape_dim = size(BFMmodel.shapePC, 2);
  msz.n_tex_dim   = size(BFMmodel.texPC,   2);
  msz.n_seg       = size(BFMmodel.segbin,  2); 
end  