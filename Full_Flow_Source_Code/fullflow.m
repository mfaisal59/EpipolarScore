function flow=fullflow(im1_ori,im2_ori,ratio,maxDisp,opt)
	im1=max(min(imresize(im1_ori(1:end-mod(size(im1_ori,1),ratio),1:end-mod(size(im1_ori,2),ratio),:),1/ratio),1),0);
	im2=max(min(imresize(im2_ori(1:end-mod(size(im2_ori,1),ratio),1:end-mod(size(im2_ori,2),ratio),:),1/ratio),1),0);
    [m,n,ch]=size(im1);
    range=-floor(maxDisp/ratio):floor(maxDisp/ratio);
	dcol=reshape(repmat(range,[numel(range) 1]),1,[]);
	drow=reshape(repmat(range',[1 numel(range)]),1,[]);
	win=1;
	label=numel(drow);		
	ma=range(end)+win;
	unary=opt.outOfRange*ones(label,m,n,'single');
	parfor i=win+1:m-win
		for j=win+1:n-win
			template=im1(i-win:i+win,j-win:j+win,:);
            if(sum(sum(abs(template(1,1,1)-template(:,:,1))))<1e-4||sum(sum(abs(template(1,1,2)-template(:,:,2))))<1e-4||sum(sum(abs(template(1,1,3)-template(:,:,3))))<1e-4)
				continue;
            end
			flag=zeros(numel(range),numel(range));
			target=im2(max(i-ma,1):min(i+ma,m),max(j-ma,1):min(j+ma,n),:);
			flag(1+max(i-ma,1)-(i-ma):end-(i+ma-min(i+ma,m)),1+max(j-ma,1)-(j-ma):end-(j+ma-min(j+ma,n)))=1;
			tmp=(1-(normxcorr2(template(:,:,1),target(:,:,1))+normxcorr2(template(:,:,2),target(:,:,2))+normxcorr2(template(:,:,3),target(:,:,3)))/3)/2;
            buf=unary(:,i,j);
            buf(find(flag))=min(reshape(tmp(1+2*win:end-2*win,1+2*win:end-2*win),1,[]),0.5);
			unary(:,i,j)=buf;
		end
	end
	unary=single(reshape(unary,label,m*n));
	a=[find([ones(m-1,n);zeros(1,n)])' find([ones(m,n-1) zeros(m,1)])'];
	b=[find([ones(m-1,n);zeros(1,n)])'+1 find([ones(m,n-1) zeros(m,1)])'+m];
	c=exp(-sqrt((im1(a)-im1(b)).^2+(im1(a+m*n)-im1(b+m*n)).^2+(im1(a+m*n*2)-im1(b+m*n*2)).^2)*opt.inverse_b);
	pairwise=single(abs(repmat(drow,[label 1])-repmat(drow,[label 1])')+abs(repmat(dcol,[label 1])-repmat(dcol,[label 1])'));
    tic
	[f,energy]=trws_l1(unary,a,b,single(c*opt.lambda),pairwise,opt.maxIter,min(max(pairwise(:)),opt.truncation));
    toc
	flow=cat(3,reshape(dcol(f),[m n]),reshape(drow(f),[m n]));
end
