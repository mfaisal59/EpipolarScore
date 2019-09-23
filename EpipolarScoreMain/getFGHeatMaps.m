function getFGHeatMaps(W, costs, imRepo)

dataset = imRepo.dataset;
imgPath = imRepo.imgPath;
outPath = imRepo.fgHeatMap;

scale = imRepo.scale;
ht = imRepo.ht;
wd = imRepo.wd;

dr = imRepo.dr;
T = length(dr);

th = prctile(costs(costs>0),95);
costs = 250*costs/th;

numsuperpixels = round((ht*wd)/400);    % 400 pixel superpixel
compactness = 20.0;

for i=1:T
    imName = sprintf('%s%s/%s', imgPath, dataset, dr(i).name);
    imgc = imresize(imread(imName), scale);
    
    [labels, ~] = snic_mex(imgc,numsuperpixels,compactness);
    spIdx = mylabel2idx(labels+1);
    
    xy = round(W(2*i-1:2*i, :));
    indc_1 = ht*(xy(1,:)-1) + xy(2,:);
    vld = ~isnan(indc_1) & xy(1,:)>0 & xy(2,:)>0 & xy(1,:)<=wd & ...
        xy(2,:)<=ht;
    
    disImg = nan(ht, wd);
    disImg(indc_1(vld)) = costs(vld);

    % find the median distance in each superpixel
%     spxlAvg = cellfun(@(x)nanmedian(disImg(x)), spIdx);
%     disImgMedian = spxlUnary2Img([ht,wd], spIdx, spxlAvg);
%     
%     figure(1)
%     subplot(2,1,1)
%     imagesc(disImg)
%     title(i)
%     axis image
%     
%     subplot(2,1,2)
%     imagesc(disImgMedian);
%     axis image
%     colormap gray
%     pause(0.01);
    
    outName = sprintf('%s%s/%s', outPath, dataset, dr(i).name);
    imwrite(uint8(disImg), outName);
end