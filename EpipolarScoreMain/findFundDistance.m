function costs = findFundDistance(Wh, Fs, frms)

T = size(Wh,1)/3;
P = size(Wh,2);

bigF = zeros(3*T);
n = size(frms,1);
for i=1:n
    f1 = frms(i,1);
    f2 = frms(i,2);
    F12 = Fs(:,:,i);
    
    bigF(3*f1-2:3*f1, 3*f2-2:3*f2) = F12';
    bigF(3*f2-2:3*f2, 3*f1-2:3*f1) = F12;
end

epipolarCost = zeros(1,P);
vstCnt = zeros(1,P);
visMsk = Wh(1:3:end, :)~=0;
for i=1:T
    Fi = bigF(:, 3*i-2:3*i);
    vldI = visMsk(i,:);
    xi = Wh(3*i-2:3*i, vldI);
    Li = Fi*xi;
   
    mskI = visMsk(:, vldI);
    mskI(i,:) = 0;
    Lnorm = sqrt(Li(1:3:end, :).^2 + Li(2:3:end, :).^2);
%     Lnorm(i,:) = 1;
    Lnorm(~mskI) = 1;
    
    Pi = size(xi,2);
    Lnh = zeros(3*T, Pi);
    Lnh(1:3:end,:) = Lnorm;
    Lnh(2:3:end,:) = Lnorm;
    Lnh(3:3:end,:) = Lnorm;
    
    Li = Li./Lnh;
    WLi = Wh(:, vldI).*Li;
    epipolarCost(vldI) = epipolarCost(vldI) + sum(abs(WLi(1:3:end,:) + WLi(2:3:end,:) + WLi(3:3:end,:)),1);
    vstCnt(vldI) = vstCnt(vldI) + sum(mskI,1);
end
vstCnt(vstCnt==0) = 1;
costs = epipolarCost./vstCnt;   % Average epipolar cost
