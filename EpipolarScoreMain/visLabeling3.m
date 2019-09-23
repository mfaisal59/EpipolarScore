function visLabeling3(W, costs, imRepo)

dataset = imRepo.dataset;
imgPath = imRepo.imgPath;
scale = imRepo.scale;
ht = imRepo.ht;
wd = imRepo.wd;
oldP = ht*wd;
% scale  = 1;
dr = imRepo.dr;
T = length(dr);

tl = 0.03;
th = 0.1;

idxb = costs<tl;
idxf = costs>th;

for i=1:T
    imName = sprintf('%s%s/%s', imgPath, dataset, dr(i).name);
    imgc = imresize(imread(imName), scale);
    
    xy = round(W(2*i-1:2*i, :));
    indc_1 = ht*(xy(1,:)-1) + xy(2,:);
    vld = ~isnan(indc_1) & indc_1>0 & indc_1<=oldP;
    
    disImg = zeros(ht, wd);
    disImg(indc_1(vld)) = costs(vld);
    
%     bgmask = false(ht, wd);
%     bgmask(indc_1(idxf&vld)) = 1;
%     imgc(:,:,1) = imgc(:,:,1) + uint8(80*bgmask);
% 
%     bgmask = false(ht, wd);
%     bgmask(indc_1(idxb&vld)) = 1;
%     imgc(:,:,2) = imgc(:,:,2) + uint8(80*bgmask);
    
    figure(1)
    subplot(1,2,1)
    imagesc(imgc)
    title(i)
    axis image
    
    subplot(1,2,2)
%     mesh(disImg);
    imagesc(disImg);
    axis image
    colormap gray
    pause;
end