function T = trifocal6pt(x1,x2,x3)
% Estimation of trifocal tensor from six points. Algorithm 20.1 p511.
% Hartley and Zisserman (2nd Ed)
% Arguments:
%          x1, x2, x3 - Three sets of corresponding 3x6 set of homogeneous
%          points.
%
% Returns:
%          T      - A cell array of size 3 or 1 with 3x3x3 trifocal 
%                   tensor inside such that 
%                   x1^i * x2^j * x3^k * eps_jpr * eps_kqs * T_i^pq = 0
%                   OR [x2]_x * SUM_i(x1^i*T_i) * [x3]_x = 0_(3x3)
%

if nargin==1
    x = x1;
    x1 = x(1:3,:);
    x2 = x(4:6,:);
    x3 = x(7:9,:);
end
% x1 = rand(2,6);
% x2 = rand(2,6);
% x3 = rand(2,6);

N = size(x1,2);
if N~=6
   error('This function requires exact six number of correspondances'); 
end
if size(x1,1)==2
    x1 = [x1; ones(1,N)];
    x2 = [x2; ones(1,N)];
    x3 = [x3; ones(1,N)];
end

% select 4 noncolliner points
[x1, x2, x3, fndsubset] = orderPoints(x1, x2, x3);
if ~fndsubset
    T = zeros(3,3,3);
    warning('Couldnt find a subset of 4 with no three collinear points');
    return;
end
E = [1 0 0; 0 1 0; 0 0 1; 1 1 1]'; % [e1, e2, e3, e4]

% should I normalize the 4 points??
l = inv(x1(:,3:5))*x1(:,6);
T1 = inv([l(1)*x1(:,3), l(2)*x1(:,4), l(3)*x1(:,5)]);

l = inv(x2(:,3:5))*x2(:,6);
T2 = inv([l(1)*x2(:,3), l(2)*x2(:,4), l(3)*x2(:,5)]);

l = inv(x3(:,3:5))*x3(:,6);
T3 = inv([l(1)*x3(:,3), l(2)*x3(:,4), l(3)*x3(:,5)]);

% E1 = T1*x1(:,3:6);
% E2 = T2*x2(:,3:6);
% E3 = T3*x3(:,3:6);

% transform the remaining two points
x1h = [T1*x1(:,1), T2*x2(:,1), T3*x3(:,1)];
x2h = [T1*x1(:,2), T2*x2(:,2), T3*x3(:,2)];

% should I normalize x1h and x2h ??
F = redcdFundamentalMatrix(x1h, x2h);

X = [1 1 1 1; 0 0 0 0; 1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]';
for i=1:length(F)
    
    tF = [F{i}(1,2) F{i}(2,1) 0; F{i}(1,3) 0 F{i}(3,1); 0 F{i}(2,3) F{i}(3,2)];
    ln = null(tF);
    if isempty(ln)
        [~,~,V] = svd(tF);
        ln = V(:,3);
    end
    rn = null(F{i}');
    if isempty(rn)
        [~,~,V] = svd(F{i}');
        rn = V(:,3);
    end
%     norm(tF*ln)
%     norm(rn'*F{i})
    
    A = [ln(2) -ln(1) 0 0; ln(3) 0 -ln(1) 0; ...
        rn(2) -rn(1) 0 (rn(1)-rn(2)); rn(3) 0 -rn(1) (rn(1)-rn(3))];
    a = null(A);
    if isempty(a)
        [~,~,V] = svd(A);
        a = V(:,4);
    end
%     norm(tF*a(1:3))
%     norm((a(4)-a(1:3))'*F{i})
    X(:,2) = a;
    
    % normalization??
    P1 = getCameraMatrix(X,x1);
    P2 = getCameraMatrix(X,x2);
    P3 = getCameraMatrix(X,x3);
    T{i} = trifocalFromP(P1, P2, P3);
    
%     [F21,F31,F32] = trifocal2FundMatrix(T{i});
%     err1 = sum(x2.*(F21*x1), 1);
%     err2 = sum(x3.*(F31*x1), 1);
%     err3 = sum(x3.*(F32*x2), 1);
    
%     verifyTrifocal(T{i}, x1, x2, x3);
end

% Find camera matrices using X and x1,x2,x3 using DLT
% Find trifocal tensor from camera matrices

function F = redcdFundamentalMatrix(x1, x2)
% Reduced Fundamental matrix, the one that corresponds to cameras
% P=[eye(3), ones(3,1)] and P' = [diag([a,b,c]), d*ones(3,1)]. 
% See Hartley, Zisserman (2nd Ed), p509.

A = ones(4, 6);
ind = [2:4, 6:8];
for i=1:3
    x1tx2t = kron(x1(:,i)', x2(:,i)');
    A(i,:) = x1tx2t(ind);
end
f = null(A)';

a(1) = f(2,1)*f(2,4)*f(2,5) + f(2,2)*f(2,3)*f(2,6);
a(2) = f(1,1)*f(2,4)*f(2,5) + f(1,4)*f(2,1)*f(2,5) + ...
    f(1,5)*f(2,1)*f(2,4) + f(1,2)*f(2,3)*f(2,6) + ...
    f(1,3)*f(2,2)*f(2,6) + f(1,6)*f(2,2)*f(2,3);
a(3) = f(1,1)*f(1,4)*f(2,5) + f(1,1)*f(1,5)*f(2,4) + ...
    f(1,4)*f(1,5)*f(2,1) + f(1,2)*f(1,3)*f(2,6) + ...
    f(1,2)*f(1,6)*f(2,3) + f(1,3)*f(1,6)*f(2,2);
a(4) = f(1,1)*f(1,4)*f(1,5) + f(1,2)*f(1,3)*f(1,6);

lamdas = roots(a);

n = 1;
for i=1:length(lamdas)
    if isreal(lamdas(i))
        fi = f(1,:) + lamdas(i)*f(2,:);
        F{n} = [0 fi(1) fi(2); fi(3) 0 fi(4); fi(5) fi(6) 0]';
        n = n + 1;
    end
end
% if length(lamdas)>1
%     F = F{1};
% end

function T = trifocalFromP(P1, P2, P3)
% get trifocal tensor from camera matrices, Eqn 15.1, page 367

H = eye(4);
H(1:3,1:3) = inv(P1(:,1:3));
H(1:3,4) = -H(1:3,1:3)*(P1(1:3,4));

P1 = P1*H;      % Canonical camera P=[I|0]
P2 = P2*H;
P3 = P3*H;

T = zeros(3,3,3);
for i=1:3
    T(:,:,i) = P2(:,i)*(P3(:,4))' - P2(:,4)*(P3(:,i))';
end

function [x1, x2, x3, fndsubset] = orderPoints(x1, x2, x3)
% order x1,x2,x3 such that in the last four points, no three should
% be collinear

N = size(x1,2);
th = 1e-6;
fndsubset = false;

iter = 1;
maxIter = 50;
A = zeros(3);
t = zeros(3);
while ~fndsubset
    ind1 = randsample(N,4);
    t(:,1) = ind1(1:3);
    t(:,2) = ind1(2:4);
    t(:,3) = ind1([4,1,2]);
    for i=1:3
        A(1,i) = det(x1(:,t(:,i)));
        A(2,i) = det(x2(:,t(:,i)));
        A(3,i) = det(x3(:,t(:,i)));
    end
    if all(all(abs(A)>th))
        fndsubset = true;
    end
    if iter>maxIter
        ind1=[];
        break;
    end
    iter = iter + 1;
end
ind2 = setdiff(1:6, ind1);
ind = [ind2, ind1'];
x1 = x1(:, ind);
x2 = x2(:, ind);
x3 = x3(:, ind);

function verifyTrifocal(T, x1, x2, x3)

N = size(x1,2);
Zs = zeros(N,1);
for i=1:N
    y1 = x1(:,i);
    y2 = x2(:,i);
    y3 = x3(:,i);
    
    smY1T = zeros(3);
    for j=1:3
        smY1T = smY1T + y1(j)*T(:,:,j);
    end
    Z = crossX(y2)*smY1T*crossX(y3);
    Zs(i) = norm(Z);
end

function xc = crossX(x)

xc = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
