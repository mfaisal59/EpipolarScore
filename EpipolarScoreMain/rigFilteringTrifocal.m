function [costs1] = rigFilteringTrifocal4(W, imRepo)

T = size(W,1)/2;
P = size(W,2);

height = imRepo.ht;
width = imRepo.wd;
oldP = height*width;
W(W==0) = nan;      % zero is a missing value

%%
minCmnpts = round(0.8*oldP);

mxnp1 = 1e4;       % RANSAC would only consider 10K points at max
mxnp2 = 2e3;
Ts = zeros(3,3,3,T);
frms = zeros(3,T);
n = 1;
disp('Estimating trifocal tensors ...');
for i=1:T-2
    pts1 = W(2*i-1:2*i, :);
    vld1 = ~isnan(pts1(1,:));
    for j=i+1%:T-1
%         if(staticFrames(j))
%             continue;
%         end
        pts2 = W(2*j-1:2*j, :);
        vld2 = ~isnan(pts2(1,:));
        for k=j+1%:T
%             if(staticFrames(k))
%                 continue;
%             end
%             fprintf(1, 'Estimating Trifocal Tensor for %d-%d-%d ...\n', i,j,k);
            pts3 = W(2*k-1:2*k, :);
            vld3 = ~isnan(pts3(1,:));
            cmn = find(vld1 & vld2 & vld3);
            
            npts = length(cmn);
            if npts<minCmnpts && k>(j+1)
                break;
            end
            
            if n>1
                sltdPts2 = intersect(sltdPts0,cmn);     % reuse previous inliers
            end
            
            if n>1 && length(sltdPts2)>mxnp2
                % I am not doing RANSAC again. It'd find trifocal from
                % inliers
                x1 = pts1(:,sltdPts2);
                x2 = pts2(:,sltdPts2);
                x3 = pts3(:,sltdPts2);
                
                if(isStationary(x1,x2,x3))
                    Ti = zeros(3,3,3);
%                 inliers = 1:size(x1,2);
                else
                    [Ti, inliers] = ransacfittrifocal2(x1, x2, x3, 0.005);
                end
            else
                if npts>mxnp1            % use a few random points
                    sltdPts = randomsample(npts, mxnp1);
                    sltdPts2 = cmn(sltdPts);
                else                % use all the points
                    sltdPts = 1:npts;
                    sltdPts2 = cmn(sltdPts);
                end
                x1 = pts1(:,sltdPts2);
                x2 = pts2(:,sltdPts2);
                x3 = pts3(:,sltdPts2);
                
                [Ti, inliers] = ransacfittrifocal2(x1, x2, x3, 0.005);
                sltdPts0 = sltdPts2(inliers);
            end
%             [Ti,~, ~,~] = ransacfittrifocal(x1,x2,x3,0.005);
            Ts(:,:,:,n) = Ti;
            frms(:,n) = [i,j,k];
            n = n + 1;
            
%             epiDist([i,j,k],cmn) = epiDist([i,j,k],cmn) + repmat(abs(d), 3,1);
%             vstCnt(cmn) = vstCnt(cmn) + 3;
        end
        if npts<minCmnpts
            break;
        end
    end
end
frms = frms(:,1:n-1);
disp('Estimating Epipolar Distances for trajectories...');
costs1 = findtrifocalDistance3(W,Ts,frms, imRepo);
