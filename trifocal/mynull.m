function a = mynull(A)

a = null(A);
if isempty(a)
    [~,~,V] = svd(A);
    a = V(:,end);
end
