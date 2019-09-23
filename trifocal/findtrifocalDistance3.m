function [costs] = findtrifocalDistance3(W,Ts,frms, imRepo)
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
vstCnt = zeros(1, P);
% epiDist = zeros(F, P);
costs = zeros(1,P);
for m=1:M
    i = frms(1,m);
    j = frms(2,m);
    k = frms(3,m);
    
    pts1 = W(2*i-1:2*i, :);
	vld1 = ~isnan(pts1(1,:));

	pts2 = W(2*j-1:2*j, :);
	vld2 = ~isnan(pts2(1,:));
        
    pts3 = W(2*k-1:2*k, :);
    vld3 = ~isnan(pts3(1,:));
    
%     fprintf(1, 'Estimating Trifocal Distance for %d-%d-%d ...\n', i,j,k);
    T = Ts(:,:,:,m);
    
    cmn = vld1 & vld2 & vld3;
    x1 = pts1(:, cmn);
    x2 = pts2(:, cmn);
    x3 = pts3(:, cmn);
    
    npts = size(x1,2);
    x1 = normalise2dpts([x1; ones(1,npts)]);
    x2 = normalise2dpts([x2; ones(1,npts)]);
    x3 = normalise2dpts([x3; ones(1,npts)]);
    
    [F21,F31,F32] = trifocal2FundMatrix(T);
    
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
    
    d = x2tFx1+x1tFtx2 + x3tFx1+x1tFtx3 + x3tFx2+x2tFtx3;
% %     epiDist([i,j,k],cmn) = epiDist([i,j,k],cmn) + repmat(d,3,1);
%     epiDist(i,cmn) = epiDist(i,cmn) + d;
    costs(cmn) = d + costs(cmn);
    vstCnt(cmn) = vstCnt(cmn) + 2;
end

vstCnt(vstCnt==0) = 1;
% costs = sum(epiDist,1)./vstCnt;
costs = costs./vstCnt;


% ht = imRepo.ht;
% wd = imRepo.wd;
% oldP = ht*wd;
% win = 5;
% runAvg = zeros(1,P);
% vlds = false(win,P);
% for i=1:F
%     xy = round(W(2*i-1:2*i, :));
%     indc_1 = ht*(xy(1,:)-1) + xy(2,:);
%     j = rem(i-1,win)+1;
%     if i>win
%         runAvg(vlds(j,:)) = runAvg(vlds(j,:)) - epiDist(i-win,vlds(j,:));
%     end
%     vlds(j,:) = ~isnan(indc_1) & indc_1>0 & indc_1<=oldP;
%     runAvg(vlds(j,:)) = runAvg(vlds(j,:)) + epiDist(i,vlds(j,:));
%     mxAvgCost = max(runAvg, [],1);
% end
% % 
% vstCnt(vstCnt==0) = 1;
% avgcost = sum(epiDist,1)./vstCnt;
% % mxcost = max(epiDist, [],1);
% % 
% % sumcost = costs/max(costs);
% % avgcost = avgcost/max(avgcost);
% % mxcost = mxcost/max(mxcost);
% % mxAvgCost = mxAvgCost/max(mxAvgCost);
% % epiDist = epiDist/max(max(epiDist));
% % % costs2 = costs./vstCnt;
% for i=1:F
%     imName = sprintf('%s%s/%s', imRepo.imgPath, imRepo.dataset, ...
%         imRepo.dr(i).name);
%     img1 = imresize(imread(imName), imRepo.scale);
%     
%     xy = round(W(2*i-1:2*i, :));
%     indc_1 = ht*(xy(1,:)-1) + xy(2,:);
%     vld = ~isnan(indc_1) & indc_1>0 & indc_1<=oldP;
%     
%     disImg1 = zeros(ht, wd);
%     disImg1(indc_1(vld)) = epiDist(i,vld);
% %     disImg1(indc_1(vld)) = sum(epiDist(:,vld), 1);
%     disImg2 = zeros(ht, wd);
% %     disImg2(indc_1(vld)) = sumcost(vld);
% %     disImg2(indc_1(vld)) = mxAvgCost(vld);
%     disImg2(indc_1(vld)) = runAvg(vld);
%     
%     disImg3 = zeros(ht, wd);
%     disImg3(indc_1(vld)) = avgcost(vld);
% 
%     figure(1)
%     subplot(2,2,1)
%     imagesc(img1);
%     title(i)
%     axis image
%     subplot(2,2,2)
%     imagesc(disImg1)
%     colormap(gray)
%     title('per-frame cost');
%     axis image
%     
%     subplot(2,2,3)
%     imagesc(disImg2)
%     title('Max Cost')
%     colormap(gray)
%     axis image
%     subplot(2,2,4)
%     imagesc(disImg3)
%     colormap(gray)
%     title('Mean Cost')
%     axis image
%     pause;
% end

