function Img = spxlUnary2Img(siz, spIdx, c)

lens = cellfun(@length, spIdx);
cArr = repelem(c, lens);
Idx = cell2mat(spIdx);

Img = zeros(siz);
Img(Idx) = cArr;