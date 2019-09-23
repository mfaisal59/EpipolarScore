function d = findEpiDistance(F, x1, x2)

npts = size(x1,2);
x1 = [x1; ones(1,npts)];
x2 = [x2; ones(1,npts)];
[x1, T1] = normalise2dpts(x1);
[x2, T2] = normalise2dpts(x2);

x2tFx1 = sum(x2.*(F*x1), 1);

Fx1 = F*x1;
Ftx2 = F'*x2;

% Evaluate distances
d =  x2tFx1.^2 ./ ...
    (Fx1(1,:).^2 + Fx1(2,:).^2 + Ftx2(1,:).^2 + Ftx2(2,:).^2);
	