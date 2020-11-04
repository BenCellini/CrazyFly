function [angle,m,pts,k] = tracktip(img, mask, rot, norm, npts, mode, arg)
%% get_vid_props: get properties of video
%
%   INPUT:
%       img   	:   image to track tip
%       rot    	:   rotation point
%       norm   	:   normalize contrast based on entire image or just ROI
%       npts 	:   # points to track
%       mode  	:   median or k-means
%       arg    	:   percentiles for median or # clusters for k-means
%       debug	:   show comparison for each frame
%
%   OUTPUT:
%       angle 	: measured angle [°]
%       pts   	: point locations
%       bw      : binarized image
%

if nargin < 7 || isempty(arg)
    if nargin < 6 || isempty(mode)
        mode = 'dist';
    end
    switch mode
        case 'dist' % average of distribution tails
            arg = [15 85]; % upper and lower 10 percentiles default
        case 'clust' % average of k-means cluster centroids
            arg = 2; % 2 clusters default
            %dthresh = 6; % angle STD limit at tip
    end
end
dthresh = 10; % angle STD limit at tip

% Convert to greyscale if needed
rot = double(rot);
if size(img,3) > 1
    img = rgb2gray(img);
end

if norm == 2 % normalize tracking ROI only
    img = imadjust(img, stretchlim(img(mask)));
elseif norm == 3 % normalize whole image
   img = imadjust(img);
end

% Binarize & get image reigon properties
% bw = imbinarize(img(:,:,1), thresh) & mask;
bw = imbinarize(img(:,:,1)) & mask;
%SE = strel('diamond', 3);
%bw = imerode(bw, SE);
p = regionprops(bw, 'Area', 'PixelList');

if ( isempty(p) || all(bw == mask, 'all') ) || (~any(bw(:)~=0) || npts <=1)
    angle = nan; % null tracking if the ROI is whited or blacked out
    warning('Empty ROI')
else
    % Only use largest object
    [~, ix] = max([p.Area]);
    x = p(ix).PixelList(:,1);
    y = p(ix).PixelList(:,2);           

    % Convert to polar coordinates and sort points by their radius
    [theta, rho] = cart2pol(x-rot(1), y-rot(2));
    theta = wrapTo360(rad2deg(theta));
    [~,rix] = sort(rho, 'descend');
    if length(rho) < npts
        npts = length(rho);
        warning(['# points tracked = ' num2str(npts)])
    end
    
    % Only keep the points with largest radii
    rix = rix(1:npts);
    theta = theta(rix);
    rho = rho(rix);
    
    % Convert back to cartesian coordinates & store points
    x = x(rix);
    y = y(rix);
    pts = [x, y];
    
    % Get angle from points
    switch mode
        case 'dist' % distribution tails
            k = nan;
            m = nan;
            angle = median([prctile(theta, arg(1)), prctile(theta, arg(2))]);
        case 'clust' % agglomerative clustering to find tip(s)
            if arg == 1 % no clustering to do if there's only 1 group
                angle = median(theta);
            else
                % Cluster data into specified # of groups
                %k = kmeans([theta rho], arg, 'Distance', 'cosine', 'MaxIter', 500);
                Y = pdist([rho theta], 'euclidean');
                Z = linkage(Y, 'average');
                k = cluster(Z, 'maxclust', arg);
                m = nan(arg,1);
                spread = nan(arg,1);
                k_n = nan(arg,1);
                for n = 1:arg % STD & median of each cluster
                    spread(n) = std(theta(k==n));
                    m(n) = median(theta(k==n));
                    k_n(n) = sum(k==n);
                end

                % Make sure clusters are grouped tighly (we found the true tip)
                % then remove extra clutsers
                spread_check = spread > dthresh; % test clusters for variance
                new_flag = false;
                p = 1;
                while any(spread_check)
                    warning('Likely more clusters than set')
                   	new_flag = true;
                    n_clust = arg + p; % add cluster

                    % Cluster data into specified # of clusters + p
                    k = cluster(Z, 'maxclust', n_clust);
                    m = nan(n_clust,1);
                    for n = 1:n_clust % STD & median of each cluster
                        spread(n) = std(theta(k==n));
                        m(n) = median(theta(k==n));
                    end
                    spread_check = spread > dthresh;

                    % Remove cluster farthest from center if extra cluster(s) was found
                    [~, sI] = sort(abs(m - 270), 'ascend'); % sort median by distance from center 
                    keepI = any(k == sI(1:2)', 2); % keep the two closest clusters
                    k = k(keepI);
                    k = findgroups(k); % start clusters from index 1 if 1st was removed
                    new_theta = theta(keepI);
                    new_rho = rho(keepI);
                    new_pts = pts(keepI,:);
                    m = nan(arg,1);
                    for n = 1:arg % median of each cluster
                        m(n) = median(theta(k==n));
                    end

                    p = p + 1;
                    if p > 5 % don't do this forever
                       warning('Clusters not found, breaking loop')
                       break
                    end
                end
                
                if new_flag
                   pts = new_pts;
                   theta = new_theta;
                   rho = new_rho;
                end
                
                % Sort clusters
                [m, sI] = sort(m, 'ascend'); % sort median by location
                new_k = k;
                new_k(k==sI(1)) = 1;
                new_k(k==sI(2)) = 2;
                k = new_k;
                
                % Mean angle
                angle = mean(m); % output angle is mean of each groups median
                
                if new_flag
                   angle = nan; 
                end

            end
    end
end