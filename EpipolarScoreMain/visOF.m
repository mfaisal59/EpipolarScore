function visOF(imRepo)

dataset = imRepo.dataset;
imgPath = imRepo.imgPath;
floPathFwd = imRepo.floPathFwd;
floPathBwd = imRepo.floPathBwd;

scale = imRepo.scale;
typ = imRepo.typ;

dr = imRepo.dr;
T = length(dr);

i = 1;
imName = sprintf('%s%s/%s', imgPath, dataset, dr(i).name);
img1 = rgb2gray(imresize(imread(imName), scale));
[height, width, ~] = size(img1);

for i=1:T-1
    imName = sprintf('%s%s/%s', imgPath, dataset, dr(i+1).name);
    img2 = rgb2gray(imresize(imread(imName), scale));

    floFwdName = sprintf('%s%s/%s.flo', floPathFwd, dataset, dr(i).name(1:end-4));
    floBwdName = sprintf('%s%s/%s.flo', floPathBwd, dataset, dr(i).name(1:end-4));
    uv = readFlowFile(floFwdName);
    uvinv = readFlowFile(floBwdName);
    
    uv = imresize(uv, scale);
    uvinv = imresize(uvinv, scale);
    
    [xc_p, yc_p] = findConsistentflow(uv, uvinv);
    
    
    ind2_1 = height*(round(xc_p)-1) + round(yc_p);
    
    indOb2_1 = ~isnan(ind2_1) & xc_p>0 & yc_p>0 & xc_p<=width & yc_p<height;
    img2_1 = 255*ones(size(img2));
    img2_1(indOb2_1) = img2(ind2_1(indOb2_1));
    
    figure(1)
    imagesc(img1)
    colormap(gray)
    figure(2)
    imagesc(img2_1)
    colormap(gray)
    
    pause(0.1);
    img1 = img2;
end
