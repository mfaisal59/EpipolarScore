function [xc_1, yc_1] = nxtFrameMapping(xp_1, yp_1, xc_p, yc_p)

[height, width] = size(xc_p);
xpl = floor(xp_1);
ypl = floor(yp_1);
xpu = ceil(xp_1);
ypu = ceil(yp_1);
xpr = round(xp_1);
ypr = round(yp_1);

xdu = xp_1 - xpl;
xdl = xpu - xp_1;
ydu = yp_1 - ypl;
ydl = ypu - yp_1;

wll = sqrt(xdl.^2 + ydl.^2);
wuu = sqrt(xdu.^2 + ydu.^2);
wlu = sqrt(xdl.^2 + ydu.^2);
wul = sqrt(xdu.^2 + ydl.^2);

% wll = xdl + ydl;
% wuu = xdu + ydu;
% wlu = xdl + ydu;
% wul = xdu + ydl;

indpll = height*(xpl-1) + ypl;
indpuu = height*(xpu-1) + ypu;
indplu = height*(xpl-1) + ypu;
indpul = height*(xpu-1) + ypl;
indpr = height*(xpr-1) + ypr;

wsm = wll+wuu+wlu+wul;      % should be = 4

vld = ~( isnan(xp_1) | isnan(yp_1) | xpl<1 | xpu>width | ...
    ypl<1 | ypu>height) ;

xc_1 = nan(size(xp_1));
yc_1 = nan(size(yp_1));

xc_1(vld) = ( wll(vld).*xc_p(indpll(vld)) + wlu(vld).*xc_p(indplu(vld)) + ...
    wul(vld).*xc_p(indpul(vld)) + wuu(vld).*xc_p(indpuu(vld)) )./wsm(vld);
ivldx = isnan(xc_1) | ~isfinite(xc_1);

yc_1(vld) = ( wll(vld).*yc_p(indpll(vld)) + wlu(vld).*yc_p(indplu(vld)) + ...
    wul(vld).*yc_p(indpul(vld)) + wuu(vld).*yc_p(indpuu(vld)) )./wsm(vld);
ivldy = isnan(yc_1) | ~isfinite(yc_1);

xc_1t = nan(size(xp_1));
vldr = ~( isnan(xpr) | isnan(ypr) | xpr<1 | xpr>width | ...
    ypr<1 | ypr>height) ;
xc_1t(vldr) = xc_p(indpr(vldr));
xc_1(ivldx) = xc_1t(ivldx);


yc_1t = nan(size(xp_1));
yc_1t(vldr) = yc_p(indpr(vldr));
yc_1(ivldy) = yc_1t(ivldy);


