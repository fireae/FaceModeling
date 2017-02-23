function [R,lambda] = linreg( X , Y , lambda )

%X = [ones(size(X,1),1) X];

%% method 1: soving linear regression using close-form solution %%%%%%%%%%%

%R = (X'*X+eye(size(X,2))*lambda)\X'*Y;

%% method 2: using SVR in liblinear %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
featdim   = size(X,2);
shapedim  = size(Y,2);

param = sprintf('-s 12 -p 0 -c %f -q', lambda);
%-s:12-- L2-regularized L2-loss support vector regression (dual)
%-p:0-- epsilon : epsilon-SVR loss function parameter-epsilon（default0.1）
%-c
%-q: quiet mode:no output information
%param = sprintf('-s 12 -p 0 -c 0.3 -q');
R_tmp = zeros( featdim, shapedim );
tic;
for o = 1 : shapedim
    disp(['Training landmarks ' num2str(o)]);
    model = train(Y(:,o),sparse(X),param);
    R_tmp(:,o) = model.w';
end
toc;

R = R_tmp;

end
