% Demo for Structured Edge Detector (please see readme.txt first).
addpath(genpath('../toolbox'));
%% set opts for training (see edgesTrain.m)
opts=edgesTrain();                % default options (good settings)
opts.modelDir='models/';          % model will be in models/forest
opts.modelFnm='modelBsds';        % model name
opts.nPos=5e5; opts.nNeg=5e5;     % decrease to speedup training
opts.useParfor=0;                 % parallelize if sufficient memory

%% train edge detector (~20m/8Gb per tree, proportional to nPos/nNeg)
tic, model=edgesTrain(opts); toc; % will load model if already trained

%% set detection parameters (can set after training)
model.opts.multiscale=0;          % for top accuracy set multiscale=1
model.opts.sharpen=2;             % for top speed set sharpen=0
model.opts.nTreesEval=4;          % for top speed set nTreesEval=1
model.opts.nThreads=4;            % max number threads for evaluation
model.opts.nms=0;                 % set to true to enable nms

%% evaluate edge detector on BSDS500 (see edgesEval.m)
if(0), edgesEval( model, 'show',1, 'name','' ); end

%% detect edge and visualize results
for id=13%[89 90 91 93 94]
	mkdir(sprintf('F:/KITTI/dataset/EpicFlow_v1.00/iphone/%04d/edges',id));
	for i=0:10000
		if(~exist(sprintf('F:/KITTI/dataset/EpicFlow_v1.00/iphone/%04d/images_sampled/%06d.png',id,i),'file'))
			continue;
		end
		i
		tic
		I = imread(sprintf('F:/KITTI/dataset/EpicFlow_v1.00/iphone/%04d/images_sampled/%06d.png',id,i));
		edges=edgesDetect(I,model);
		%save(sprintf('F:/KITTI/dataset/edges/%02d/image_2/%06d.mat',id,i),'edges');
		fid=fopen(sprintf('F:/KITTI/dataset/EpicFlow_v1.00/iphone/%04d/edges/%06d.bin',id,i),'wb');
		fwrite(fid,transpose(edges),'single');
		fclose(fid);		
		toc
	end
end
