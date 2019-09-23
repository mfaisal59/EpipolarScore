% set up:
% add path for trifocal folder and compile snic_mex.cpp
% addpath('../OF+MS/mycodes/external/trifocal');
%
% create a root dirctory containing forward and backward flow and 
% the images. These folders further contain subfolders for each 
% of the sequences. Also create the folder, 'fgHeatMap'. This would
% contain the heat maps for your seqeunces.

imRepo.root = '/home/cvlab/FAISAL/DAVIS_Dataset/DAVIS-trainval-480p/';
imRepo.floPathBwd = [imRepo.root, 'fullFlow/davisbackward/'];
imRepo.floPathFwd = [imRepo.root, 'fullFlow/davisforward/'];
imRepo.imgPath = [imRepo.root, 'JPEGImages/480p/'];
imRepo.fgHeatMap = [imRepo.root, 'epipolarScore/480p/'];

imRepo.scale = [480, 854];     % resize the image to make it tractable
% This is the same scale we used when estimating optical flow
imRepo.typ = 'jpg';

drp = dir(imRepo.imgPath);      % A list of all the seqeunces
for k=3:length(drp) % check k should start with 3 or 4 (depends on OS)
    imRepo.dataset = drp(k).name;   % pick a dataset
%     imRepo.dataset = 'bear';
    disp(imRepo.dataset);
    mkdir([imRepo.fgHeatMap,imRepo.dataset]);
    
    cmd = sprintf('%s%s/*.%s', imRepo.imgPath, imRepo.dataset, imRepo.typ);
    imRepo.dr = dir(cmd);       % all the images in the dataset
    visOF(imRepo)
end
