function bool = isStationary(x1, x2, x3)

N = size(x1,2);
minStcPts = round(0.5*N);
th = 1e-3;

d12 = abs(x1-x2);
d13 = abs(x1-x3);
d23 = abs(x2-x3);

flgs = (d12(1,:)<th & d12(2,:)<th) | (d13(1,:)<th & d13(2,:)<th) | ...
    (d23(1,:)<th & d23(2,:)<th);
staticPts = sum(flgs);

bool = staticPts>minStcPts;