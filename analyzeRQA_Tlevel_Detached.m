function OUT = analyzeRQA_Tlevel_Detached(RQA, Y, cx, cy, opts)

if nargin < 5
    opts = struct();
end

if ~isfield(opts,'minGapPx'), opts.minGapPx = 50; end
if ~isfield(opts,'nearCentroidEndPx'), opts.nearCentroidEndPx = 80; end
if ~isfield(opts,'dilateSeedPx'), opts.dilateSeedPx = 1; end
if ~isfield(opts,'minKeepArea'), opts.minKeepArea = 500; end

RemoveSeed = false(size(Y));
detachedRay = false(size(RQA.theta));
detachedSeg = cell(size(RQA.theta));
kList = [];
removeStartR = [];
removeEndR = [];

for k = 1:numel(RQA.theta)

    s = RQA.segStartR{k};
    e = RQA.segEndR{k};

    if numel(s) < 2
        detachedSeg{k} = [];
        continue
    end

    gap = s(2:end) - e(1:end-1);

    remove_j = false(size(s));

    for j = 1:numel(s)

        % Core T-level detached rule:
        % segment is very close to centroid
        isNearCentroid = e(j) < opts.nearCentroidEndPx;

        % has large separation from next muscle-like segment
        hasLargeGap = false;
        if j < numel(s)
            hasLargeGap = gap(j) >= opts.minGapPx;
        end

        if isNearCentroid && hasLargeGap
            remove_j(j) = true;
            kList(end+1) = k;
            removeStartR(end+1) = s(j);
            removeEndR(end+1) = e(j);
        end
    end

    detachedSeg{k} = find(remove_j);

    if any(remove_j)
        detachedRay(k) = true;

        t = RQA.theta(k);

        for jj = find(remove_j(:))'

            rr_remove = s(jj):e(jj);

            x = round(cx + rr_remove .* cos(t));
            y = round(cy + rr_remove .* sin(t));

            valid = x >= 1 & x <= size(Y,2) & ...
                    y >= 1 & y <= size(Y,1);

            idx = sub2ind(size(Y), y(valid), x(valid));
            RemoveSeed(idx) = true;
        end
    end
end
% Extend detected angular support to nearby nYS==1 rays
%{
if ~isfield(opts,'extendSingleWindow'), opts.extendSingleWindow = 5; end

coreK = unique(kList);
for kk0 = coreK(:)'

    neigh = kk0-opts.extendSingleWindow : kk0+opts.extendSingleWindow;
    neigh = mod(neigh-1, numel(RQA.theta)) + 1;

    for kk = neigh(:)'
        if RQA.nYSegments(kk) == 1 && ~ismember(kk,kList)

            s1 = RQA.segStartR{kk};
            e1 = RQA.segEndR{kk};

            if ~isempty(s1) && e1(1) < opts.nearCentroidEndPx
                kList(end+1) = kk;
                removeStartR(end+1) = 1;
                removeEndR(end+1) = min(e1(1), opts.nearCentroidEndPx);
            end
        end
    end
end
%} 
% Slight dilation to ensure seed touches the leakage component

if opts.dilateSeedPx > 0
    RemoveSeed = imdilate(RemoveSeed, strel('disk', opts.dilateSeedPx));
end

if ~isempty(kList)

    kList0 = kList(:);
    rStart0 = removeStartR(:);
    rEnd0   = removeEndR(:);

    % handle circular wrap: 355-360 and 1-10 should be one cluster
    kList2 = kList0;
    kList2(kList2 <= 10) = kList2(kList2 <= 10) + 360;

    [ks, ord] = sort(kList2);
    rs = rStart0(ord);
    re = rEnd0(ord);

    gapK = diff(ks);
    %breakIdx = [0; find(gapK > 3); numel(ks)];
    breakIdx = [0; find(gapK > opts.extendSingleWindow+1); numel(ks)];
    RemoveRegion = false(size(Y));

    for c = 1:numel(breakIdx)-1

        id = breakIdx(c)+1 : breakIdx(c+1);

        kk = mod(ks(id)-1, 360) + 1;
        thetaList = RQA.theta(kk);

        rIn  = ones(size(re(id))) * 1; %zeros(size(re(id)));  %0;  rs(id);
        rOut = re(id);

        xIn  = cx + rIn  .* cos(thetaList(:));
        yIn  = cy + rIn  .* sin(thetaList(:));

        xOut = cx + rOut .* cos(thetaList(:));
        yOut = cy + rOut .* sin(thetaList(:));

        xPoly = [xIn; flipud(xOut)];
        yPoly = [yIn; flipud(yOut)];

        RemoveRegion = RemoveRegion | ...
            poly2mask(xPoly, yPoly, size(Y,1), size(Y,2));
    end

else
    RemoveRegion = false(size(Y));
end

RemoveRegion = RemoveRegion & Y;

% Marker-controlled leakage removal
%{
Ytmp = Y & ~RemoveRegion;
CC = bwconncomp(Ytmp);

RemoveComponent = false(size(Y));

for c = 1:CC.NumObjects
    
    pix = CC.PixelIdxList{c};

    if any(RemoveRegion(pix))
        RemoveComponent(pix) = true;
    end

end
%}
Ytmp = Y & ~RemoveRegion;   % ตัด sector ก่อน

CC = bwconncomp(Ytmp);
stats = regionprops(CC,'Area');

Yrev = false(size(Y));

for c = 1:CC.NumObjects
    if stats(c).Area >= opts.minKeepArea
        Yrev(CC.PixelIdxList{c}) = true;
    end
end

%Yrev = Ytmp & ~RemoveComponent;

%% ...OUT structure ...
OUT = struct();
OUT.detachedRay = detachedRay;
OUT.detachedTheta = RQA.thetaRef(detachedRay);
OUT.detachedSeg = detachedSeg;
OUT.RemoveSeed = RemoveSeed;
OUT.RemoveRegion = RemoveRegion;
OUT.Yrev = Yrev;
OUT.kList = kList;
OUT.removeStartR = removeStartR;
OUT.removeEndR = removeEndR;
OUT.RemoveRegion = RemoveRegion;
OUT.Yrev = Yrev;
end