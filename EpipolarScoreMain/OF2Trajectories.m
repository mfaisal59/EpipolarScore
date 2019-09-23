function [W] = OF2Trajectories4(imRepo, drctn)

% Due to occlusion, new pixels are assigned new labels gives in newPxlLbls 
% pxLblsInv takes back from new labels to traditional pixel labels (indRf).
% W gives the mapping from current frame to first (or first obsreved pixel)
% Columns in W are the 2D projection of trajectories of points in 3D
% map1st2Crnt is opposite to W and gives image coordinates that maps first 
% image to the current frame. The pixels not present in the first frame are 
% mapped to the frame where they were observed the very first time. 

dataset = imRepo.dataset;
imgPath = imRepo.imgPath;
scale = imRepo.scale;
if nargin<2
    drctn = 1;
end
T = length(imRepo.dr);
ht = imRepo.ht;
wd = imRepo.wd;
P = ht*wd;

[xr,yr] = meshgrid(1:wd, 1:ht);
indRf = ht*(xr(:)'-1) + yr(:)';        % rest (reference) coordinates

Xc_p = nan(ht, wd, T-1);
Yc_p = nan(ht, wd, T-1);

if(drctn>0)
    floPathFwd = imRepo.floPathFwd;
    floPathBwd = imRepo.floPathBwd;
    dr = imRepo.dr;
    
    for i=2:T
        floFwdName = sprintf('%s%s/%s.flo', floPathFwd, dataset, dr(i-1).name(1:end-4));
        floBwdName = sprintf('%s%s/%s.flo', floPathBwd, dataset, dr(i-1).name(1:end-4));
        
        [~, ~, Xc_p(:,:,i-1), Yc_p(:,:,i-1)] = getConsistentFlow(...
            floFwdName, floBwdName,scale);
    end
else
    floPathBwd = imRepo.floPathFwd;
    floPathFwd = imRepo.floPathBwd;
    dr = flipud(imRepo.dr);
    
    for i=2:T
        floFwdName = sprintf('%s%s/%s.flo', floPathFwd, dataset, dr(i).name(1:end-4));
        floBwdName = sprintf('%s%s/%s.flo', floPathBwd, dataset, dr(i).name(1:end-4));
        
        [~, ~, Xc_p(:,:,i-1), Yc_p(:,:,i-1)] = getConsistentFlow(...
            floFwdName, floBwdName, scale);
    end
end

%%
oldP = P;
W = nan(2*T, P);
W(1:2, 1:oldP) = [xr(:), yr(:)]';
xc_1 = xr(:)';
yc_1 = yr(:)';

[newPxlLbls,~] = meshgrid(1:P, 1:T);
pxLblsInv = zeros(T, P);
pxLblsInv(1,:) = indRf;
for i=2:T        
    xp_1 = xc_1;
    yp_1 = yc_1;
    xc_p = Xc_p(:,:,i-1);
    yc_p = Yc_p(:,:,i-1);
    
    rxp_1 = round(xp_1);
    ryp_1 = round(yp_1);
    indp_1 = ht*(rxp_1-1) + ryp_1;
    ocPxp = isnan(indp_1) | rxp_1<1 | ryp_1<1 | rxp_1>wd | ryp_1>ht;
    xc_1(:) = nan;
    yc_1(:) = nan;
    
    [xc_1(~ocPxp), yc_1(~ocPxp)] = nxtFrameMapping(xp_1(~ocPxp), yp_1(~ocPxp), xc_p, yc_p);
    
    xc = round(xc_1);
    yc = round(yc_1);
    indc_1 = ht*(xc-1) + round(yc);
    vldc_1 = ~isnan(indc_1) & xc>0 & yc>0 & xc<=wd & yc<=ht;
%     [~,ia] = unique(indc_1, 'stable');
%     uflg = false(size(vldc_1));
%     uflg(ia) = 1;
%     vldc_1 = vldc_1 & uflg;
%     xc_1(~uflg) = 0;
%     yc_1(~uflg) = 0;
%     vldc_1 = ~isnan(indc_1) & indc_1>0 & indc_1<=oldP;
    lstPxl = setdiff(indRf, indc_1(vldc_1(:))');
    
    npx = length(lstPxl);
    newPxlLbls(i, lstPxl) = P+1:P+npx;
    P = P + npx;
    pxLblsInv(i, newPxlLbls(i,:)) = indRf;
    xc_1(newPxlLbls(i, lstPxl)) = xr(lstPxl);
    
    yc_1(newPxlLbls(i, lstPxl)) = yr(lstPxl);
    pc = length(xc_1);
    W(2*i-1:2*i, 1:pc) = [xc_1; yc_1];
end
W(W==0) = nan;

function [uv, uvinv, xc_p, yc_p] = getConsistentFlow(floFwdname, floBwdname, scale)
% function [uv, uvinv] = getGroundFlow(floName, occName, scale)

uv = readFlowFile(floFwdname);
uvinv = readFlowFile(floBwdname);
uv = imresize(uv, scale);
uvinv = imresize(uvinv, scale);
[xc_p, yc_p] = findConsistentflow(uv, uvinv);