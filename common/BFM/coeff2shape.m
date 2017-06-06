function shape = coeff2shape(CrTensor, wID,wEP)
%shape [3 N] 
%Re_shape  = coef2object(para, model.shapeMU, model.shapePC, model.shapeEV );
scale2mm = 100;
shape = scale2mm*double(ttm(CrTensor,{wID,wEP},[2,3]));
shape = reshape(shape,3,[]);

end