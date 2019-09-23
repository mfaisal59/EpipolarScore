function [spIdx2, ulbls, idx, cnts] = mylabel2idx(lbls)
% Find indices for each lable
% input labels should be >0, else labels<1 would be ignored
% output Idx is a cell array of size=max(labels), consisting of
% the pixel indices in input labels. The 2nd output is the unique 
% labels (sorted) present in input labels
% 

lbls(isnan(lbls)) = 0;
[ulbls, cnts, idx] = myunique(lbls(:));

% idx = (1:numel(lbls));
% [srtdLbls, ind] = sort(lbls(:));
% idx = idx(ind);
% 
% [ulbls, ia] = unique(srtdLbls, 'rows', 'stable');
% ia = [ia; size(srtdLbls,1)+1];
% cnts = ia(2:end) - ia(1:end-1);

spIdx = mat2cell(idx, 1, cnts);
vld = ulbls>0;
ulbls = ulbls(vld);
spIdx = spIdx(vld);
spIdx2 = cell(1,ulbls(end));
% spIdx2 = cell(1,srtdLbls(end));
spIdx2(ulbls) = spIdx;

% N = max(max(lbls));
% for i=1:N
%    spIdx2{i} = find(lbls==i)';
% end
% diff = cellfun(@minus, spIdx, spIdx2, 'uni', 0);
% diff2 = cellfun(@(x)power(x,2), diff, 'uni', 0);
% err = cellfun(@sum, diff2);