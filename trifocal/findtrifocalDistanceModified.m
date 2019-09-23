function [h] = findtrifocalDistanceModified(W,Ts,frms, imRepo, vstCnt)
% Find symetric epipolar distance of a point set (x1 or x2 or x3) 
% from the epiplolar lines im the two other views. costs is the 
% mean of all three (rather six) distances
%
% W consisst of the traajectories, Ts is 4D array of trifocal 
% tensors, the corresponding frame #s are given in frms

F = size(W,1)/2;
P = size(W,2);

if nargin >2
    M = size(frms,2);
else
    M = 1;
    frms = [1 2 3]';
end

WW = zeros(3*M,P);
FF = zeros(3*M,3*M);

dif=0;
rows = 0;
cols = 4;
for f=1:F
    pts1 = W(2*f-1:2*f, :);
    vld1 = ~isnan(pts1(1,:));
    x1 = pts1(:,vld1);
    npts = size(x1,2);
    temp = normalise2dpts([x1; ones(1,npts)]);
    WW(2*f-1+dif:2*f+dif, vld1) = temp(1:2,:);
    WW(3*f,vld1) = temp(3,:);
%     WW(2*m-1+dif:2*m+dif, :) = W(2*m-1:2*m, :);
%     WW(3*m,:) = 1;
    dif = dif + 1;
end

for m = 1 : M
    i = frms(1,m);
    j = frms(2,m);
    k = frms(3,m);
    T = Ts(:,:,:,m);
    [F12,F13,F23] = trifocal2FundMatrix(T);
    
    FF(i*3-2:i*3, j*3-2:j*3) = F12;
    FF(j*3-2:j*3, i*3-2:i*3) = F12;
    
    FF(i*3-2:i*3, (j+1)*3-2:(j+1)*3) = F13;
    FF((j+1)*3-2:(j+1)*3, i*3-2:i*3) = F13;
    
    FF((i+1)*3-2:(i+1)*3, (j+1)*3-2:(j+1)*3) = F23;
    FF((i+1)*3-2:(i+1)*3, j*3-2:j*3) = F23;  
end
vstCnt(vstCnt==0) = 1;
h = zeros(1,P);
h = abs(sum(WW.*(FF*WW),1));
h = h./vstCnt;
% vstCnt = zeros(1, P);
% % epiDist = zeros(F, P);
% costs = zeros(1,P);
% for m=1:M-2
%     i = frms(1,m);
%     j = frms(2,m);
%     k = frms(3,m);
%     
%     pts1 = W(2*i-1:2*i, :);
% 	vld1 = ~isnan(pts1(1,:));
% 
% 	pts2 = W(2*j-1:2*j, :);
% 	vld2 = ~isnan(pts2(1,:));
%         
%     pts3 = W(2*k-1:2*k, :);
%     vld3 = ~isnan(pts3(1,:));
%     
% %     fprintf(1, 'Estimating Trifocal Distance for %d-%d-%d ...\n', i,j,k);
%     T = Ts(:,:,:,m);
%     
%     cmn = vld1 & vld2 & vld3;
%     x1 = pts1(:, cmn);
%     x2 = pts2(:, cmn);
%     x3 = pts3(:, cmn);
%     
%     npts = size(x1,2);
%     x1 = normalise2dpts([x1; ones(1,npts)]);
%     x2 = normalise2dpts([x2; ones(1,npts)]);
%     x3 = normalise2dpts([x3; ones(1,npts)]);
%     
%     [F21,F31,F32] = trifocal2FundMatrix(T);
%     
%     l21 = F21*x1;
%     l12 = F21'*x2;
%     
%     l31 = F31*x1;
%     l13 = F31'*x3;
%     
%     l32 = F32*x2;
%     l23 = F32'*x3;
%     
%     l21 = l21./repmat(sqrt(l21(1,:).^2 + l21(2,:).^2), 3,1);
%     l12 = l12./repmat(sqrt(l12(1,:).^2 + l12(2,:).^2), 3,1);
%     
%     l31 = l31./repmat(sqrt(l31(1,:).^2 + l31(2,:).^2), 3,1);
%     l13 = l13./repmat(sqrt(l13(1,:).^2 + l13(2,:).^2), 3,1);
%     
%     l32 = l32./repmat(sqrt(l32(1,:).^2 + l32(2,:).^2), 3,1);
%     l23 = l23./repmat(sqrt(l23(1,:).^2 + l23(2,:).^2), 3,1);
%     
%     x2tFx1 = abs(sum(x2.*l21, 1));
%     x1tFtx2 = abs(sum(x1.*l12, 1));
%     
%     x3tFx1 = abs(sum(x3.*l31, 1));
%     x1tFtx3 = abs(sum(x1.*l13, 1));
%     
%     x3tFx2 = abs(sum(x3.*l32, 1));
%     x2tFtx3 = abs(sum(x2.*l23, 1));
%     
%     d = x2tFx1+x1tFtx2 + x3tFx1+x1tFtx3 + x3tFx2+x2tFtx3;
% % %     epiDist([i,j,k],cmn) = epiDist([i,j,k],cmn) + repmat(d,3,1);
% %     epiDist(i,cmn) = epiDist(i,cmn) + d;
%     costs(cmn) = d + costs(cmn);
%     vstCnt(cmn) = vstCnt(cmn) + 2;
% end
% 
% vstCnt(vstCnt==0) = 1;
% % costs = sum(epiDist,1)./vstCnt;
% costs = costs./vstCnt;