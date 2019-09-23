function visW(W, imRepo, Ig)

scale = imRepo.scale;
dataset = imRepo.dataset;
imgPath = imRepo.imgPath;
dr = imRepo.dr;
if size(Ig,1)==3
    Ig = 0.21*Ig(1,:) + 0.72*Ig(2,:) + 0.07*Ig(3,:);
end
i=1;
imName = sprintf('%s%s/%s', imgPath, dataset, dr(i).name);
img1 = rgb2gray(imresize(imread(imName), scale));

height = size(img1,1);
width = size(img1,2);
[xr,yr] = meshgrid(1:width, 1:height);

oldP = height*width;
T = size(W,1)/2;
P = size(W,2);
W(W==0) = nan;

% Iw = nan(T,P);
% for i=1:T
%     imName = sprintf('%s%s/%s%04d.%s', imgPath, dataset, pre, i, typ);
%     imgc = rgb2gray(imresize(imread(imName), scale));
%     
%     xy = round(W(2*i-1:2*i, :));
%     indc_1 = height*(xy(1,:)-1) + xy(2,:);
%     vld = ~isnan(indc_1) & indc_1>0 & indc_1<=oldP;
%     Iw(i,vld) = imgc(indc_1(vld));
% end
% visPts = sum(~isnan(Iw),1);
% Iw(isnan(Iw))=0;
% Ig = sum(Iw,1)./visPts;

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

MSGID = 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId';
warning('off', MSGID);
% Ig(Ig==0)= 255;
for i=1:T
    imName = sprintf('%s%s/%s', imgPath, dataset, dr(i).name);
    imgc = rgb2gray(imresize(imread(imName), scale));
    
    xy = round(W(2*i-1:2*i, :));
    indc_1 = height*(xy(1,:)-1) + xy(2,:);
    vld = ~isnan(indc_1) & indc_1>0 & xy(1,:)<=width & xy(2,:)<=height;
    
    imgc_1 = 255*zeros(height, width);
    imgc_1(indc_1(vld)) = Ig(vld);
    
%     xy = W(2*i-1:2*i, :);
%     vld = ~isnan(xy(1,:));
%     Itc_1 = scatteredInterpolant(xy(1,vld)', xy(2,vld)', Ig(vld)', 'natural','nearest');
%     imgc_1 = Itc_1(xr, yr);
   
    figure(1)
    imagesc(imgc)
    colormap(gray)
    axis image
    figure(3)
    imagesc(imgc_1)
    colormap(gray)
    axis image
    pause;
end
warning('on', MSGID);