% function Fs = rigFiltering

% addpath('../floReadWrite');
% load OFTrajecories_v1_scale0.5b.mat

W(W==0) = nan;
dataset = 'alley_2';
imgPath = '../../../Sintel/MPI-Sintel-training_images/training/clean/';

scale = 1;
typ = 'png';
pre = 'frame_';
cmd = sprintf('%s%s/*.%s', imgPath, dataset, typ);
dr = dir(cmd);

f = 1;
% F = size(U,1);
imName = sprintf('%s%s/%s%04d.%s', imgPath, dataset, pre, f, typ);
img = imresize(imread(imName), scale);
[height,width,~] = size(img);

T = length(dr);

imName = sprintf('%s%s/%s%04d.%s', imgPath, dataset, pre, 1, typ);
img1 = imresize(imread(imName), scale);

[height, width, ~] = size(img1);
oldP = height*width;

pts1 = W(1:2, :);
vld1 = ~isnan(pts1(1,:));
pts2 = W(3:4, :);
vld2 = ~isnan(pts2(1,:));

for i=3:T
    imName = sprintf('%s%s/%s%04d.%s', imgPath, dataset, pre, i, typ);
    img = imresize(imread(imName), scale);

    disp(sprintf('Trifocal Tensor for %d-%d-%d ...\n', i-2, i-1, i));
    
    pts3 = W(2*i-1:2*i, :);
    vld3 = ~isnan(pts3(1,:));
    
    cmn = vld1 & vld2 & vld3;
    
%     [~,~, ~,~, d] = ransacfittrifocal(pts1(:, cmn), pts2(:, cmn), ...
%         pts3(:,cmn), 0.01);
    
%     [~, d2] = ransacfittrifocal2(pts1(:, cmn), pts2(:, cmn), ...
%         pts3(:,cmn), 0.01);

        [Ti] = ransacfittrifocal2(pts1(:, cmn), pts2(:, cmn), ...
        pts3(:,cmn), 0.01);
    [d] = findtrifocalDistance2(W(2*i-5:2*i, cmn), Ti);
    [d2] = findtrifocalDistance3(W(2*i-5:2*i, cmn), Ti);
    
    indc_1 = height*(round(pts3(1,:))-1) + round(pts3(2,:));
    vld = ~isnan(indc_1) & indc_1>0 & indc_1<=oldP;
    
    err = nan(height, width);
    err(indc_1(cmn)) = d;
    img1 = img;
    mask = err<2;
    img1(:,:,2) = img1(:,:,2) + uint8(40*mask);
    
    figure(1)
    subplot(2,2,1)
    imagesc(img1);
    axis image
    subplot(2,2,2)
    imagesc(err)
    axis image
    
    err2 = nan(height, width);
    err2(indc_1(cmn)) = d2;
    img2 = img;
    mask = err<2;
    img2(:,:,2) = img2(:,:,2) + uint8(40*mask);
    
    subplot(2,2,3)
    imagesc(img2);
    axis image
    subplot(2,2,4)
    imagesc(err2)
    axis image
    
%     figure(2)
%     for j=1:6
%         err1 = nan(height, width);
%         err1(indc_1(cmn)) = d2{j};
%         
%         subplot(3,2,j)
%         imagesc(err1)
%         axis image
%     end
    
    
    pause(0.01);
    
    pts1 = pts2;
    vld1 = vld2;
    pts2 = pts3;
    vld2 = vld3;
%     img1 = img2;
end