function [costs] = findtrifocalDistance2(W,Ts,frms)
% Find symetric epipolar distance of a point set (x1 or x2 or x3) 
% from the epiplolar lines im the two other views. costs is the 
% mean of all three (rather six) distances
%
% W consisst of the traajectories, Ts is 4D array of trifocal 
% tensors, the corresponding frame #s are given in frms

T = size(W,1)/2;
P = size(W,2);
if nargin >2
    M = size(frms,2);
else
    M = 1;
    frms = [1 2 3]';
end
vstCnt = zeros(1, P);
epiDist = zeros(T, P);

for m=1:M
    i = frms(1,m);
    j = frms(2,m);
    k = frms(3,m);
    
    fprintf(1, 'Estimating Trifocal Distance for %d-%d-%d ...\n', i,j,k);
    
    T = Ts(:,:,:,m);
    [F21,F31,F32] = trifocal2FundMatrix(T);
    
    if m==1 || frms(2,m-1)~=frms(2,m)
        pts1 = W(2*i-1:2*i, :);
        vld1 = ~isnan(pts1(1,:));
        pts2 = W(2*j-1:2*j, :);
        vld2 = ~isnan(pts2(1,:));
    end
    pts3 = W(2*k-1:2*k, :);
    vld3 = ~isnan(pts3(1,:));
    
    if m==1 || frms(2,m-1)~=frms(2,m)
        cmn = vld1 & vld2 & vld3;
        cmn2 = cmn;
        x1 = pts1(:, cmn);
        x2 = pts2(:, cmn);
        x3 = pts3(:, cmn);
        
        npts = size(x1,2);
        x1 = [x1; ones(1,npts)];
        x2 = [x2; ones(1,npts)];
        x3 = [x3; ones(1,npts)];
        
        l21 = F21*x1;
        l12 = F21'*x2;
        l31 = F31*x1;
        l13 = F31'*x3;
        l32 = F32*x2;
        l23 = F32'*x3;
        
        l21 = l21./repmat(sqrt(l21(1,:).^2 + l21(2,:).^2), 3,1);
        l12 = l12./repmat(sqrt(l12(1,:).^2 + l12(2,:).^2), 3,1);
        l31 = l31./repmat(sqrt(l31(1,:).^2 + l31(2,:).^2), 3,1);
        l13 = l13./repmat(sqrt(l13(1,:).^2 + l13(2,:).^2), 3,1);
        l32 = l32./repmat(sqrt(l32(1,:).^2 + l32(2,:).^2), 3,1);
        l23 = l23./repmat(sqrt(l23(1,:).^2 + l23(2,:).^2), 3,1);
        
        x2tFx1 = abs(sum(x2.*l21, 1));
        x1tFtx2 = abs(sum(x1.*l12, 1));
        x3tFx1 = abs(sum(x3.*l31, 1));
        x1tFtx3 = abs(sum(x1.*l13, 1));
        x3tFx2 = abs(sum(x3.*l32, 1));
        x2tFtx3 = abs(sum(x2.*l23, 1));
    else
        cmn2 = cmn & vld3;
        
        x3 = pts3(:, cmn);      % it also includes the lost pixels
        x3 = [x3; ones(1,npts)];
        lst = find(isnan(x3(1,:)));
        x3(:,lst) = 0;
        
        l31 = F31*x1;
        l13 = F31'*x3;
        l32 = F32*x2;
        l23 = F32'*x3;
        
        l31 = l31./repmat(sqrt(l31(1,:).^2 + l31(2,:).^2), 3,1);
        l13 = l13./repmat(sqrt(l13(1,:).^2 + l13(2,:).^2), 3,1);
        l32 = l32./repmat(sqrt(l32(1,:).^2 + l32(2,:).^2), 3,1);
        l23 = l23./repmat(sqrt(l23(1,:).^2 + l23(2,:).^2), 3,1);
        
        x2tFx1(lst) = 0;
        x1tFtx2(lst) = 0;
        x3tFx1 = abs(sum(x3.*l31, 1));
        x1tFtx3 = abs(sum(x1.*l13, 1));
        x3tFx2 = abs(sum(x3.*l32, 1));
        x2tFtx3 = abs(sum(x2.*l23, 1));
    end
    
    d = x2tFx1+x1tFtx2 + x3tFx1+x1tFtx3 + x3tFx2+x2tFtx3;
    epiDist([i,j,k],cmn) = epiDist([i,j,k],cmn) + repmat(d,3,1);
    vstCnt(cmn2) = vstCnt(cmn2) + 6;
end
vstCnt(vstCnt==0) = 1;
costs = sum(epiDist,1)./vstCnt;