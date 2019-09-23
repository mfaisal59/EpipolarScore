function [C, cnts, idx] = myunique(A)
% Give unique elements and their counts. Idx contains the indices 
% of these elements in A. unique(A(idx)) = C
% Input should be numeric 1D array

idx = 1:length(A);
[srtdA, ind] = sort(A);
idx = idx(ind);
[C, ia] = unique(srtdA, 'rows', 'stable');
ia = [ia; size(srtdA,1)+1];
cnts = ia(2:end) - ia(1:end-1);