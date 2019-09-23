% RANSACFITTRIFOCAL - fits trifocal tensor using RANSAC
%
% Usage:   [T, e2, e3, inliers, d] = ransacfittrifocal(x1, x2, x3, t)
%
% Arguments:
%          x1  - 2xN or 3xN set of homogeneous points.  If the data is
%                2xN it is assumed the homogeneous scale factor is 1.
%          x2,x3  - 2xN or 3xN set of homogeneous points such that 
%               x1<->x2<->x3
%          t   - The distance threshold between data point and the model
%                used to decide whether a point is an inlier or not. 
%                Note that point coordinates are normalised to that their
%                mean distance from the origin is sqrt(2).  The value of
%                t should be set relative to this, say in the range 
%                0.01 - 0.1  
%
% Note that it is assumed that the matching of x1, x2 & x3 are putative and
% it is expected that a percentage of matches will be wrong.
%
% Returns:
%          T       - The 3x3x3 trifocal tensor such that 
%                   x1^i * x2^j * x3^k * eps_jpr * eps_kqs * T_i^pq = 0
%                   OR [x2]_x * SUM_i(x1^i*T_i) * [x3]_x = 0_(3x3)
%          e2, e3 - The epipoles in image 2 and image 3, corresponding to  
%                   the first image camera center
%          inliers - An array of indices of the elements of x1, x2 that were
%                    the inliers for the best model.
%          d        - the trifocal distance of the points x1,x2 & x3
%


function [T, inliers] = ransacfittrifocal2(x1, x2, x3, th, inliers)

    if ~all(size(x1)==size(x2)) || ~all(size(x1)==size(x3))
        error('Data sets x1, x2 & x3 must have the same dimension');
    end
    
    [rows,npts] = size(x1);
    if rows~=2 && rows~=3
        error('x1 and x2 must have 2 or 3 rows');
    end
    
    if rows == 2    % Pad data with homogeneous scale factor of 1
        x1 = [x1; ones(1,npts)];
        x2 = [x2; ones(1,npts)];
        x3 = [x3; ones(1,npts)];
    end
    
    % Normalise each set of points so that the origin is at centroid and
    % mean distance from origin is sqrt(2).  normalise2dpts also ensures the
    % scale parameter is 1.
    
%     meandist = mean(sqrt(x1(1,:).^2 + x1(2,:).^2));
%     scale = sqrt(2)/meandist;
%     th = scale*th;
    
    [x1, H1] = normalise2dpts(x1);
    [x2, H2] = normalise2dpts(x2);
    [x3, H3] = normalise2dpts(x3);
    
    s = 6;  % Number of points needed to fit a trifocal tensor for the 
            % linear algorithm

    if nargin==4        
        fittingfn = @trifocal6pt;
        distfn    = @tridist;
        degenfn   = @isdegenerate;
        % x1 and x2 are 'stacked' to create a 6xN array for ransac
        [T, inliers] = ransac([x1; x2; x3], fittingfn, distfn, degenfn, s, th);
    end
    % Now do a final least squares fit on the data points considered to
    % be inliers.
    
    [T] = triftensor(x1(:, inliers), x2(:, inliers), x3(:, inliers));
    
%     x = [x1; x2; x3];
%     [~,~,d] = tridist({T}, x, th);
%     d = d';
    % Denormalise
%     Y = zeros(3,3,3);
%     for ii = 1:3,
%         Y(:,:,ii) = inv(H2) * T(:,:,ii) * (inv(H3))';
%     end
%     T_dnrm = zeros(3,3,3);
%     for ii = 1:3,
%         T_dnrm(:,:,ii) = H1(1,ii)*Y(:,:,1) + H1(2,ii)*Y(:,:,2) ...
%                            + H1(3,ii)*Y(:,:,3);
%     end
%     T = T_dnrm;
    
function [bestInliers, bestT, d] = tridist(T, x, th)
    
    x1 = x(1:3,:);    % Extract x1, x2, & x3 from x
    x2 = x(4:6,:);
    x3 = x(7:9,:);
    npts = size(x,2);
    
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
    
    ninliers = 0;     % Number of inliers
    bestInliers = [];
    bestT = T{1};
    for i=1:length(T)
        Ti = permute( T{i} , [2 1 3] );

        t = Ti(:);
        b = A*t;
        d = abs(b(1:npts)) + abs(b(npts+1:2*npts)) + ...
            abs(b(2*npts+1:3*npts)) + abs(b(3*npts+1:4*npts));
        inliers =  find(abs(d) < th);     % Indices of inlying points
        if length(inliers) > ninliers
            ninliers = length(inliers);
            bestT = Ti;
            bestInliers = inliers;
        end
    end
%----------------------------------------------------------------------
% (Degenerate!) function to determine if a set of matched points will
% result in a degeneracy in the calculation of a trifocal tensor as 
% needed by RANSAC.  This function assumes this cannot happen...
     
function r = isdegenerate(x)
    r = 0;    
    
