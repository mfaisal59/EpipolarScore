function removeOcclusion(forward,backward,ratio,filename)
	[m,n,ch]=size(forward);
	cols=repmat(1:n,[m 1]);
	rows=repmat((1:m)',[1 n]);
	pos=cat(3,cols,rows);
	feature1=reshape(cat(3,pos,pos+forward),[],4);
	feature2=reshape(cat(3,pos+backward,pos),[],4);
	[idx,dis]=knnsearch(feature2,feature1);
	occ=exp(-reshape(dis,m,n).^2/2^2);
	threshold=0.8;
	output=fopen(filename,'w');
	for i=1:m
		for j=1:n
			if(occ(i,j)>threshold&&j+forward(i,j,1)>=1&&j+forward(i,j,1)<=n&&i+forward(i,j,2)>=1&&i+forward(i,j,2)<=m)
				fprintf(output,'%g %g %g %g\n',ratio*j-ratio/2-0.5,ratio*i-ratio/2-0.5,ratio*(j+forward(i,j,1))-ratio/2-0.5,ratio*(i+forward(i,j,2))-ratio/2-0.5);                    
			end
		end
	end
	fclose(output);
end
