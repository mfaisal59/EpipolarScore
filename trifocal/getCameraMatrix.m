function P = getCameraMatrix(X,x)
% DLT Algorithm for Camera matrix estimation from six or more points 
% correspondance (image and the world). p178, Hartley & Zisserman(2nd Ed)

% for i=1:size(X,2)
%    Xh = [X(:,i)', 1];
%    z = [0 0 0 0];
%    A(2*i-1, :) = [Xh, z, -x(1,i)*Xh];
%    A(2*i, :) = [z, Xh, -x(2,i)*Xh];
% end

N = size(X,2);
if size(X,1)==3
    X = [X; ones(1,N)];
end
if size(x,1)==2
    x = [x; ones(1,N)];
end

A = zeros(2*N, 12);
z = zeros(1,4);
for i=1:N
   A(2*i-1, :) = [z, -x(3,i)*X(:,i)', x(2,i)*X(:,i)'];
   A(2*i, :) = [x(3,i)*X(:,i)', z, -x(1,i)*X(:,i)'];
end
n = null(A);
if isempty(n)
    [~,~,V] = svd(A);
    n = V(:, end);
end
P = reshape(n, 4, 3)';