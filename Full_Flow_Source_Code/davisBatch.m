addpath(genpath('external'));%load external libraries

imgPath = '/home/cvlab/FAISAL/DAVIS_Dataset/DAVIS-trainval-480p/JPEGImages/480p/';
typ = 'jpg';
scale = [480, 854]; %adjust scale for fast compuation of OF
floPath1 = 'forward';
floPath2 = 'backward';

drp = dir(imgPath);
for k=3:length(drp)
    dataset = drp(k).name;

    mkdir(sprintf('%s/%s', floPath1, dataset));
    mkdir(sprintf('debug/forward/%s', dataset));
    mkdir(sprintf('%s/%s', floPath2, dataset));
    mkdir(sprintf('debug/backward/%s', dataset));

    cmd = sprintf('%s%s/*.%s', imgPath, dataset, typ);
    dr = dir(cmd);

    % addpath(genpath('external'));   %load external libraries
    ratio=3;        %downsample ratio
    maxDisp=100;    %smaller displacement for fast processing
    %maxDisp=242;   %maximal displacement used in the evaluation
    opt=struct('setting','sintel','outOfRange',0.22,'lambda',0.021,'truncation',1e8,'maxIter',3,'inverse_b',3);%Sintel paramters

    for i=1:length(dr)-1
        disp(sprintf('Finding flow for %d-%d...\n', i, i+1));

        name1 = sprintf('%s%s/%s', imgPath, dataset, dr(i+1).name);
        name2 = sprintf('%s%s/%s', imgPath, dataset, dr(i).name);    
        im1=im2double(imresize(imread(name1),scale));
        im2=im2double(imresize(imread(name2),scale));

        forward=fullflow(im1,im2,ratio,maxDisp,opt);        %compute the forward flow
        backward=fullflow(im2,im1,ratio,maxDisp,opt);       %compute the backward flow

        mtchFile1 = sprintf('debug/forward/%s/%s.txt', dataset, dr(i).name(1:end-4));
        removeOcclusion(forward,backward,ratio, mtchFile1);  %save the matches to match.txt after forward-backward consistency checking
        flow=runEpicflow(im1,im2, mtchFile1, opt.setting);   %apply Epicflow interpolation on the matches
        %figure(1)
        %imshow(flowToColor(flow));                          %show the flow visualization
        %pause(0.001)
        
        outName = sprintf('%s/%s/%s.flo', floPath1, dataset, dr(i).name(1:end-4));
        writeFlowFile(flow, outName);

        mtchFile2 = sprintf('debug/backward/%s/%s.txt', dataset, dr(i).name(1:end-4));
        removeOcclusion(backward,forward,ratio, mtchFile2);  %save the matches to match.txt after forward-backward consistency checking
        flow=runEpicflow(im2,im1, mtchFile2, opt.setting);   %apply Epicflow interpolation on the matches
        %figure(2)
        %imshow(flowToColor(flow));                          %show the flow visualization
        %pause(0.001)
        
        outName = sprintf('%s/%s/%s.flo', floPath2, dataset, dr(i).name(1:end-4));
        writeFlowFile(flow, outName);
    end    
end