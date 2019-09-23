% Demo for Structured Edge Detector (please see readme.txt first).
addpath(genpath('../toolbox'));
%% set opts for test (see edgesTrain.m)
opts=edgesTrain();                % default options (good settings)
opts.modelDir='models/';          % model will be in models/forest
opts.modelFnm='modelBsds';        % model name
opts.nPos=5e5; opts.nNeg=5e5;     % decrease to speedup test
opts.useParfor=0;                 % parallelize if sufficient memory

%% train edge detector (~20m/8Gb per tree, proportional to nPos/nNeg)
tic, model=edgesTrain(opts); toc; % will load model if already trained

%% set detection parameters (can set after test)
model.opts.multiscale=0;          % for top accuracy set multiscale=1
model.opts.sharpen=2;             % for top speed set sharpen=0
model.opts.nTreesEval=4;          % for top speed set nTreesEval=1
model.opts.nThreads=4;            % max number threads for evaluation
model.opts.nms=0;                 % set to true to enable nms

%% evaluate edge detector on BSDS500 (see edgesEval.m)
if(0), edgesEval( model, 'show',1, 'name','' ); end

%% detect edge and visualize results
for i=0:199
	if(~exist(sprintf('../../../data/kitti_2012/training/colored_0/better_%06d_10.png',i),'file'))
		continue;
	end
	I = imread(sprintf('../../../data/kitti_2012/training/colored_0/better_%06d_10.png',i));
	edges=edgesDetect(I,model);
	fid=fopen(sprintf('../../../data/kitti_2012/training/colored_0/better_%06d_10.bin',i),'wb');
	fwrite(fid,transpose(edges),'single');
	fclose(fid);
end

