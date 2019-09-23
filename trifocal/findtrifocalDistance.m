function costs = findtrifocalDistance(W, Ts, frms)
% forms a linear relation of [x2]_x * SUM_i(x1^i*T_i) * 
% [x3]_x = 0_(3x3) and find the value of this trilinear term
% W consisst of the traajectories, Ts is 4D array of trifocal 
% tensors, the corresponding frame #s are given in frms

T = size(W,1)/2;
P = size(W,2);
M = size(frms,2);

vstCnt = zeros(1, P);
epiDist = zeros(T, P);

for m=1:M
    i = frms(1,m);
    j = frms(2,m);
    k = frms(3,m);
    
    % else the 1st and the 2nd frame are the same as the previous
    if m==1 || frms(2,m-1)~=frms(2,m)
        pts1 = W(2*i-1:2*i, :);
        vld1 = ~isnan(pts1(1,:));

        pts2 = W(2*j-1:2*j, :);
        vld2 = ~isnan(pts2(1,:));
    end
    pts3 = W(2*k-1:2*k, :);
    vld3 = ~isnan(pts3(1,:));
    
    fprintf(1, 'Estimating Trifocal Distance for %d-%d-%d ...\n', i,j,k);
    
    
    T = Ts(:,:,:,m);
    if m==1 || frms(2,m-1)~=frms(2,m)
        cmn = vld1 & vld2 & vld3;
        
        meandist = mean(sqrt(pts1(1,cmn).^2 + pts1(2,cmn).^2));
        scale = sqrt(2)/meandist;
%         scale =1;
        
        x1 = scale*pts1(:, cmn);
        x2 = scale*pts2(:, cmn);
        x3 = scale*pts3(:, cmn);
        cmn2 = cmn;
        
        npts = size(x1,2);
        x1 = [x1; ones(1,npts)];
        x2 = [x2; ones(1,npts)];
        x3 = [x3; ones(1,npts)];
        
        A = zeros(4*npts,27);
        n = 0;
        for I = 1:2,
            for L = 1:2,
                n = n+1;
                r = (n-1)*npts+1 : n*npts;
                c1 = 3*(I-1)+L; % = [1  2  4  5]
                c2 = 3*(I-1)+3; % = [3  3  6  6]
                c3 = L+6;       % = [7  8  7  8]
                c4 = 9;         % = [9  9  9  9]
                c = [0 9 18];
                
                A(r, c1+c) =  x1';
                A(r, c2+c) = -x1' .* ( x3(L,:)' * ones(1,3) );
                A(r, c3+c) = -x1' .* ( x2(I,:)' * ones(1,3) );
                A(r, c4+c) =  x1' .* ( x2(I,:)'.*x3(L,:)' * ones(1,3) );
            end
        end
    else
        cmn2 = cmn & vld3;
        
        x3 = scale*pts3(:, cmn);      % it also includes the lost pixels
        x3 = [x3; ones(1,npts)];
        lst = find(isnan(x3(1,:)));
        x3(:,lst) = 0;          % 
        
        n = 0;
        for I = 1:2,
            for L = 1:2,
                n = n+1;
                r = (n-1)*npts+1 : n*npts;
                c2 = 3*(I-1)+3; % = [3  3  6  6]
                c4 = 9;         % = [9  9  9  9]
                c = [0 9 18];
                
                A(r, c2+c) = -x1' .* ( x3(L,:)' * ones(1,3) );
                A(r, c4+c) =  x1' .* ( x2(I,:)'.*x3(L,:)' * ones(1,3) );
                
                lr = (n-1)*npts + lst;
                A(lr,:) = 0;
            end
        end
    end
    H = (1/scale)*eye(3);
    H(3,3) = 1;
    T = permute( T , [2 1 3] );
    for ii = 1:3,
        Y(:,:,ii) = inv(H) * T(:,:,ii) * (inv(H))';
    end
    
    for ii = 1:3,
        T(:,:,ii) = H(1,ii)*Y(:,:,1) + H(2,ii)*Y(:,:,2) ...
                           + H(3,ii)*Y(:,:,3);
    end
    
    t = T(:);
    b = A*t;
    d = (abs(b(1:npts)) + abs(b(npts+1:2*npts)) + ...
        abs(b(2*npts+1:3*npts)) + abs(b(3*npts+1:4*npts)))';
    
    epiDist([i,j,k],cmn) = epiDist([i,j,k],cmn) + repmat(abs(d),3,1);
    vstCnt(cmn2) = vstCnt(cmn2) + 6;
end

vstCnt(vstCnt==0) = 1;
costs = sum(epiDist,1)./vstCnt;