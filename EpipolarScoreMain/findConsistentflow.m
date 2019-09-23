function [xc_p, yc_p] = findConsistentflow(uv, uvinv)
[height, width, ~] = size(uv);
oldP = height*width;
[xr,yr] = meshgrid(1:width, 1:height);

xc_p = xr + uv(:,:,1);
yc_p = yr + uv(:,:,2);
xp_c = xr + uvinv(:,:,1);
yp_c = yr + uvinv(:,:,2);

indc_p = height*(round(xc_p)-1) + round(yc_p);
indp_c = height*(round(xp_c)-1) + round(yp_c);
% [xc_p2, yc_p2] = getinverseFlowMapping(xp_c, yp_c);
% indc_p2 = height*(round(xc_p2)-1) + round(yc_p2);

indc_p2 = nan(size(indc_p));
[tmp, I] = sort(indp_c(:));
vld = tmp>0 & tmp<=oldP & ~isnan(tmp);
indc_p2(tmp(vld)) = I(vld);

mtchd = indc_p == indc_p2;
% mtchd = abs(indp_c-indp_c1)<5;
xc_p(~mtchd) = nan;
yc_p(~mtchd) = nan;