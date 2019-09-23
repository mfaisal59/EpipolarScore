% set up:
% add path for trifocal folder and compile snic_mex.cpp
addpath('../trifocal');
%
% create a root dirctory containing forward and backward flow and 
% the images. These folders further contain subfolders for each 
% of the sequences. Also create the folder, 'fgHeatMap'. This would
% contain the heat maps for your seqeunces.
clear;
clc;
imRepo.root = '/home/cvlab/FAISAL/DAVIS_Dataset/DAVIS-trainval-480p/';
imRepo.floPathBwd = [imRepo.root, 'fullFlow/davisbackward/'];
imRepo.floPathFwd = [imRepo.root, 'fullFlow/davisforward/'];
imRepo.imgPath = [imRepo.root, 'JPEGImages/480p/'];
imRepo.fgHeatMap = [imRepo.root, 'epipolarScore/480p/'];

imRepo.scale = [480, 854];      % resize the image to make it tractable
% This is the same scale we used when estimating optical flow
imRepo.typ = 'jpg';

drp = dir(imRepo.imgPath);      % A list of all the seqeunces
for k=3:length(drp) % check k should start with 3 or 4 (depends on OS)
    
    imRepo.dataset = drp(k).name;   % pick a dataset
    disp(imRepo.dataset);
    mkdir([imRepo.fgHeatMap,imRepo.dataset]);
    
    cmd = sprintf('%s%s/*.%s', imRepo.imgPath, imRepo.dataset, imRepo.typ);
    
    imRepo.dr = dir(cmd);       % all the images in the dataset
        
    i=1;
    imName = sprintf('%s%s/%s', imRepo.imgPath, imRepo.dataset, imRepo.dr(i).name);
    img1 = rgb2gray(imresize(imread(imName),imRepo.scale));
    
    imRepo.ht = size(img1,1);
    imRepo.wd = size(img1,2);
    
    disp('Converting Optical flow to trajectories ...');
    [W] = OF2Trajectories(imRepo);

    %get trifocal tensor based epipolar cost
    [costs] = rigFilteringTrifocal(W, imRepo);

    disp('Converting trajectory epipolar distances to epipolar score...');
    % go inside and uncomment, if you wanna visulaize/use superpixel 
    % based median/mean filtering
    getFGHeatMaps(W, costs, imRepo);    % Heatmaps in fgHeatMap folder

    disp([imRepo.dataset, ' ... done.'])
end