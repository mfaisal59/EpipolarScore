clear;
clc;
imRepo.root = './/DAVIS_Dataset/DAVIS-trainval-480p/';
imRepo.floPathBwd = [imRepo.root, 'fullFlow/davisbackward/'];
imRepo.floPathFwd = [imRepo.root, 'fullFlow/davisforward/'];
imRepo.imgPath = [imRepo.root, 'JPEGImages/'];
imRepo.fgHeatMap = [imRepo.root, 'OpticalFlow/'];

imRepo.scale = [480, 854];     % resize the image to make it tractable
% This is the same scale we used when estimating optical flow
imRepo.typ = 'jpg';

load minmax.mat;
drp = dir(imRepo.imgPath);      % A list of all the seqeunces
for k=3:length(drp) % check k should start with 3 or 4 (depends on OS)
    
    imRepo.dataset = drp(k).name;   % pick a dataset
    dataset = imRepo.dataset;
    disp(imRepo.dataset);
    cmd = sprintf('%s%s/*.%s', imRepo.imgPath, imRepo.dataset, imRepo.typ);
    imRepo.dr = dir(cmd);       % all the images in the dataset
    dr = imRepo.dr;
    T = length(dr);
    
    saveDir = sprintf('%s%s', imRepo.fgHeatMap, imRepo.dataset);
    mkdir(saveDir)
    for i=2:T
        floFwdName = sprintf('%s%s/%s.flo', imRepo.floPathFwd, dataset, dr(i-1).name(1:end-4));
        floBwdName = sprintf('%s%s/%s.flo', imRepo.floPathBwd, dataset, dr(i-1).name(1:end-4));
        
        uv = readFlowFile(floFwdName);
        uvinv = readFlowFile(floBwdName);
        
        uv = imresize(uv, imRepo.scale);
        uvinv = imresize(uvinv, imRepo.scale);
        
        u = uv(:,:,1);
        v = uv(:,:,2);
        
        u = (u-minU)/(maxU-minU);
        v = (v-minV)/(maxV-minV);
        
        u = u * 255;
        v = v * 255;
        
        fprintf('%s->%s max u = %f\n',dataset,dr(i-1).name, max(u(:)));
        fprintf('%s->%s max v = %f\n',dataset,dr(i-1).name, max(v(:)));
        
        saveNameX = sprintf('%s/flow_x_%s', saveDir, dr(i-1).name);
        saveNameY = sprintf('%s/flow_y_%s', saveDir, dr(i-1).name);
        imwrite(uint8(u), saveNameX);
        imwrite(uint8(v), saveNameY); 
    end
end


