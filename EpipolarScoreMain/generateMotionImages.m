clear
clc

datasetPath = './DAVIS_Dataset/EpipolarScore/480p';
folderPaths = dir(datasetPath);

scale = [480, 854];     % resize the image to make it tractable
% This is the same scale we used when estimating optical flow

for i = 3 : size(folderPaths,1)
    folderPath = strcat(datasetPath,'/',folderPaths(i).name);
    imagesPaths = dir(folderPath);
    outputFolderPath = strrep(datasetPath,'EpipolarScore', 'motionImages');
    mkdir(strcat(outputFolderPath,'/',folderPaths(i).name));
    
    for j = 3 : size(imagesPaths, 1)-1
        epImg = strcat(datasetPath,'/',folderPaths(i).name,'/',imagesPaths(j).name);
        heatmapImg = imresize(imread(epImg), scale);
        
        ofImgPath = strrep(datasetPath,'EpipolarScore','OpticalFlowMinMax');
        xDispImgPath = strcat(ofImgPath,'/',folderPaths(i).name,'/flow_x_',imagesPaths(j).name);
        yDispImgPath = strcat(ofImgPath,'/',folderPaths(i).name,'/flow_y_',imagesPaths(j).name);
        
        xDispImg = imresize(imread(xDispImgPath),scale);
        yDispImg = imresize(imread(yDispImgPath),scale);
        
        motionImg = zeros(size(yDispImg,1),size(yDispImg,2),3);
        motionImg(:,:,1) = yDispImg;
        motionImg(:,:,2) = xDispImg;
        motionImg(:,:,3) = heatmapImg;
            
        outputImgPath = strcat(outputFolderPath,'/',folderPaths(i).name,'/',imagesPaths(j).name);
        
        imwrite(uint8(motionImg),outputImgPath);
    end
    
    epImg = strcat(datasetPath,'/',folderPaths(i).name,'/',imagesPaths(j+1).name);
    heatmapImg = imresize(imread(epImg),[480,854]);

    ofImgPath = strrep(datasetPath,'fgHeatMapNew','OpticalFlowMinMax');
    xDispImgPath = strcat(ofImgPath,'/',folderPaths(i).name,'/flow_x_',imagesPaths(j).name);
    yDispImgPath = strcat(ofImgPath,'/',folderPaths(i).name,'/flow_y_',imagesPaths(j).name);

    xDispImg = imresize(imread(xDispImgPath),scale);
    yDispImg = imresize(imread(yDispImgPath),scale);

    motionImg = zeros(size(yDispImg,1),size(yDispImg,2),3);
    motionImg(:,:,1) = yDispImg;
    motionImg(:,:,2) = xDispImg;
    motionImg(:,:,3) = heatmapImg;
    
    outputImgPath = strcat(outputFolderPath,'/',folderPaths(i).name,'/',imagesPaths(j+1).name);

    imwrite(uint8(motionImg),outputImgPath);
end