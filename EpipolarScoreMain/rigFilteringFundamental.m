function [costs] = rigFilteringFundamental(W)

W(isnan(W)) = 0;        % Replace nans with 0
% Don't run RANSAC again and again if you already have a big set of inliers
T = size(W,1)/2;
P = size(W,2);
% numP = sum(W(1,:)==0);      % # of pixels per frame
nP = sum(W(1,:)>0);
minCmnPts = 0.1*nP;

Wh = zeros(3*T,P);
for i=1:T
    vld = W(2*i,:)>0;
    x = W(2*i-1:2*i, vld);
    
    np = size(x,2);
    xh = [x; ones(1,np)];
    
    xh = normalise2dpts(xh);
    Wh(3*i-2:3*i, vld) = xh;
end

th = 0.001;
Fs = zeros(3,3,T);
frms = zeros(T,2);
n = 1;

disp('Estimating Fundamental matrices...');
for i=1:T-1
    pts1 = Wh(3*i-2:3*i, :);
    vld1 = (pts1(1,:)~=0);
    for j=i+1:T %j=i+1:min(T, i+win)
        pts2 = Wh(3*j-2:3*j, :);
        vld2 = (pts2(1,:)~=0);
        %fprintf(1, 'Estimating Fundamental Matrix for %d-%d ...\n', i,j);
        
        cmn = find(vld1 & vld2);
        npts = length(cmn);
        if npts<minCmnPts
            break;
        end
        x1 = pts1(:,cmn);
        x2 = pts2(:,cmn);
        
        [Fi, inliers] = ransacfitfundmatrix(x1, x2, th);
        
        %%  
        Fs(:,:,n) = Fi;
        frms(n,:) = [i,j];
        n = n + 1;
    end
end
frms = frms(1:n-1,:);
disp('Estimating Epipolar Costs...');
costs = findFundDistance(Wh, Fs, frms);
