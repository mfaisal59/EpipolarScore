function [F21,F31,F32] = trifocal2FundMatrix(T)

% Estimate three pair-wise fundamental matrices given a trifocal 
% tensor. Algorirthm 15.1, page375, Hartley & Zisserman (2nd Ed) 

% [e2,e3] = e_from_T(T);
[e2,e3] = getepipoles(T);

e2 = e2/norm(e2);
e3 = e3/norm(e3);

e2x = crossX(e2);
e3x = crossX(e3);

F21 = [e2x*T(:,:,1)*e3, e2x*T(:,:,2)*e3, e2x*T(:,:,3)*e3];
F31 = [e3x*T(:,:,1)'*e2, e3x*T(:,:,2)'*e2, e3x*T(:,:,3)'*e2];

P2 = [[T(:,:,1)*e3, T(:,:,2)*e3, T(:,:,3)*e3], e2];
et = e3*e3'-eye(3);
P3 = [[et*T(:,:,1)'*e2, et*T(:,:,2)'*e2, et*T(:,:,3)'*e2], e3];

F32 = FfromPs(P2, P3);

function [e2,e3] = getepipoles(T)
% Estimate epipoles w.r.t the first frame, given a trifocal tensor
% Hartley and Zisserman p375 (2nd Ed) 

U = zeros(3);
V = zeros(3);
for i = 1:3
    % Extract epipolar lines from the column of V corresponding to
    % smallest singular value.
    U(:,i) = mynull(T(:,:,i)');
    V(:,i) = mynull(T(:,:,i));
end

% Extract epipole from the column of V corresponding to
% smallest singular value.
e2 = mynull(U');
e3 = mynull(V');

function xc = crossX(x)

xc = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];

function F = FfromPs(P1, P2)

% Estimate fundamental matrix given the pair of camera matrices.
% Eqn 9.1, page 244.  Hartley & Zisserman (2nd Ed) 

% H = eye(4);
% H(1:3,1:3) = inv(P1(:,1:3));
% H(1:3,4) = -H(1:3,1:3)*(P1(1:3,4));
% 
% P1 = P1*H;      % Canonical camera P=[I|0]
% P2 = P2*H;

C = null(P1);
e2 = P2*C;
F = crossX(e2)*P2*pinv(P1);
