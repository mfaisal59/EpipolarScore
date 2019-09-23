function Ig = getTrajColor(W, imRepo)

scale = imRepo.scale;
dataset = imRepo.dataset;
imgPath = imRepo.imgPath;
dr = imRepo.dr;
i=1;
imName = sprintf('%s%s/%s', imgPath, dataset, dr(i).name);
img1 = rgb2gray(imresize(imread(imName), scale));

height = size(img1,1);
width = size(img1,2);
% [xr,yr] = meshgrid(1:width, 1:height);
% oldP = height*width;

T = size(W,1)/2;
P = size(W,2);
W(W==0) = nan;

Iw = nan(T,P);
for i=1:T
    imName = sprintf('%s%s/%s', imRepo.imgPath, imRepo.dataset, imRepo.dr(i).name);
    imgc = rgb2gray(imresize(imread(imName), scale));
    
    xy = round(W(2*i-1:2*i, :));
    indc_1 = height*(xy(1,:)-1) + xy(2,:);
    vld = ~isnan(indc_1) & indc_1>0 & xy(1,:)<=width & xy(2,:)<=height;
%     vld = ~isnan(indc_1) & indc_1>0 & indc_1<=oldP;
    Iw(i,vld) = imgc(indc_1(vld));
end
visPts = sum(~isnan(Iw),1);
Iw(isnan(Iw))=0;
Ig = sum(Iw,1)./visPts;

% Ig = nan(1,P);
% cnt = zeros(1,P);
% for i=1:T
%     imName = sprintf('%s%s/%s%04d.%s', imgPath, dataset, pre, i, typ);
%     imgc = rgb2gray(imresize(imread(imName), scale));
%     
%     xy = round(W(2*i-1:2*i, :));
%     indc_1 = height*(xy(1,:)-1) + xy(2,:);
%     vld = ~isnan(indc_1) & indc_1>0 & xy(1,:)<=width & xy(2,:)<=height;
%     cnt(vld) = cnt(vld)+1;
%     vld2 = cnt<4 & vld;
%     Ig(vld2) = imgc(indc_1(vld2));
% end
