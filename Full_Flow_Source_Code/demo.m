addpath(genpath('external'));%load external libraries
ratio=3;%downsample ratio
maxDisp=100;%smaller displacement for fast processing
%maxDisp=242;%maximal displacement used in the evaluation
opt=struct('setting','sintel','outOfRange',0.22,'lambda',0.021,'truncation',1e8,'maxIter',3,'inverse_b',3);%Sintel paramters
%opt=struct('setting','kitti','outOfRange',0.22,'lambda',0.021,'truncation',1e8,'maxIter',3,'inverse_b',0);%Kitti paramters
im1=im2double(imread('input/a.png'));%read image 1
im2=im2double(imread('input/b.png'));%read image 2
forward=fullflow(im1,im2,ratio,maxDisp,opt);%compute the forward flow
backward=fullflow(im2,im1,ratio,maxDisp,opt);%compute the backward flow
mkdir('debug');
removeOcclusion(forward,backward,ratio,'debug/match.txt');%save the matches to match.txt after forward-backward consistency checking
flow=runEpicflow(im1,im2,'debug/match.txt',opt.setting);%apply Epicflow interpolation on the matches
imshow(flowToColor(flow));%show the flow visualization
