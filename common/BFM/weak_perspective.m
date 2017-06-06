function [R,t,s] = weak_perspective(x, X)
% reference: Optimum Fiducials under Weak Perspective Projection
% by Alfred M.Bruckstein

%获得世界坐标系到相机坐标系的弱仿射变换R和t; s为尺度
%即论文中三个矩阵:A(s),R(R),T(t)完成3D->2D的仿射变换
%s还表示尺度：2d和3d的尺度分布不同[-1,1]或者[0,width]。
%x:[2 N] X:[3 N]
% center X
mX = mean(X, 2);
X = X - repmat(mX, 1, size(X, 2));

% center x
mx = mean(x, 2);
x = x - repmat(mx, 1, size(x, 2));

% least squares
n = size(x, 2);
A = zeros(n * 2, 8);
for i = 1 : n
    A(2 * i - 1, :) = [X(:, i)', 1,  0, 0, 0, 0];
    A(2 * i, :) = [0, 0, 0, 0, X(:, i)', 1];
 end
P = A \ x(:);
P = reshape(P, 4, 2)';
lambda = P(1, 1 : 3);
gamma = P(2, 1 : 3);

% orthonormalization
d1 = sqrt((lambda * lambda') * (gamma * gamma') - (lambda * gamma') * (lambda * gamma'));
d2 = (lambda * lambda') * (gamma * gamma') + norm(lambda) * norm(gamma) * d1 - (lambda * gamma') * (lambda * gamma');
R = zeros(3);
R(1, :) = ((norm(lambda) + norm(gamma)) / (2 * norm(lambda)) + norm(gamma) * (lambda * gamma') * (lambda * gamma') / (2 * norm(lambda) * d2)) * lambda - (lambda * gamma') * gamma / (2 * d1);
R(2, :) = ((norm(lambda) + norm(gamma)) / (2 * norm(gamma)) + norm(lambda) * (lambda * gamma') * (lambda * gamma') / (2 * norm(gamma) * d2)) * gamma - (lambda * gamma') * lambda / (2 * d1);
s = norm(R(1, :)); %alpha in the paper
R(1, :) = R(1, :) / s;
R(2, :) = R(2, :) / s;
R(3, :) = cross(R(1, :), R(2, :));
t = [mx(1); mx(2); 0] - s * R * mX;
end