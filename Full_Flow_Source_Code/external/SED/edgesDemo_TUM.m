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
folder=dir('F:/TUM');
for id=3:numel(folder)
	if(~isequal(folder(id).name(1:4),'rgbd'))
		continue;
	end
	lst=dir(['F:/TUM/' folder(id).name '/rgb/*.png']);
	mkdir(['F:/TUM/edges/' folder(id).name]);
	for i=1:numel(lst)
		tic
		I=imread(['F:/TUM/' folder(id).name '/rgb/' lst(i).name]);
		edges=edgesDetect(I,model);
		fid=fopen(['F:/TUM/edges/' folder(id).name '/' lst(i).name '.bin'],'wb');
		fwrite(fid,transpose(edges),'single');
		fclose(fid);
		toc
	end
end
