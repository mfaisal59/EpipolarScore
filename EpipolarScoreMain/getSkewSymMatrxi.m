function F = getSkewSymMatrxi(a)

F = zeros(3);
F(1,2) = -a(3);
F(2,1) = a(3);

F(1,3) = a(2);
F(3,1) = -a(2);

F(2,3) = -a(1);
F(3,2) = a(1);
