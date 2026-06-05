function RQA = radialConsistencyScore(Mask_Fa, Y, opts)

if nargin < 3
    opts = struct();
end

if ~isfield(opts,'nAngles'), opts.nAngles = 180; end
if ~isfield(opts,'smoothWin'), opts.smoothWin = 9; end
if ~isfield(opts,'minSegmentPx'), opts.minSegmentPx = 3; end

Mask_Fa = logical(Mask_Fa);
Y = logical(Y);

S = regionprops(Mask_Fa,'Centroid');
C = S(1).Centroid;
cx = C(1);
cy = C(2);

[H,W] = size(Mask_Fa);

theta = linspace(0, 2*pi, opts.nAngles+1);
theta(end) = [];

bodyR = nan(size(theta));
yInner = nan(size(theta));
yOuter = nan(size(theta));
yThickness = zeros(size(theta));
nYSegments = zeros(size(theta));
maxSegLen = zeros(size(theta));
meanSegLen = zeros(size(theta));
smallSegCount = zeros(size(theta));
segmentDensity = zeros(size(theta));

maxR = hypot(H,W);

for k = 1:numel(theta)

    t = theta(k);

    rr = 1:maxR;
    xs = round(cx + rr*cos(t));
    ys = round(cy + rr*sin(t));

    valid = xs>=1 & xs<=W & ys>=1 & ys<=H;
    xs = xs(valid);
    ys = ys(valid);
    rr = rr(valid);

    idx = sub2ind([H,W],ys,xs);

    bodyLine = Mask_Fa(idx);
    yLine = Y(idx);

    % body boundary = last point inside body
    insideIdx = find(bodyLine);
    if isempty(insideIdx)
        continue;
    end

    bodyR(k) = rr(insideIdx(end));

    % restrict to inside body
    yLine(~bodyLine) = 0;

    % find Y transitions
    dY = diff([0; yLine(:); 0]);
    starts = find(dY == 1);
    ends   = find(dY == -1) - 1;

    segLen = ends - starts + 1;
    keep = segLen >= opts.minSegmentPx;

    starts = starts(keep);
    ends = ends(keep);
    segLen = segLen(keep);

    nYSegments(k) = numel(starts);
    if ~isempty(segLen)
        maxSegLen(k) = max(segLen);
        meanSegLen(k) = mean(segLen);
        smallSegCount(k) = sum(segLen < 8);
        segmentDensity(k) = sum(segLen) / max(bodyR(k),1);
   
        yInner(k) = rr(starts(1));
        yOuter(k) = rr(ends(end));
        yThickness(k) = sum(segLen);
    end
end

% Normalize radius by body radius
rInnerN = yInner ./ bodyR;
rOuterN = yOuter ./ bodyR;
thickN  = yThickness ./ bodyR;

% Smooth curves
rInnerFit = smoothdata(rInnerN,'movmedian',opts.smoothWin,'omitnan');
rOuterFit = smoothdata(rOuterN,'movmedian',opts.smoothWin,'omitnan');
thickFit  = smoothdata(thickN,'movmedian',opts.smoothWin,'omitnan');
densityFit = smoothdata(thickN,'movmedian',...
    opts.smoothWin,'omitnan');

densityResidual = thickN - densityFit;

RQA = struct();
RQA.theta = theta;
RQA.bodyR = bodyR;
RQA.rInnerN = rInnerN;
RQA.rOuterN = rOuterN;
RQA.thickN = thickN;
RQA.rInnerFit = rInnerFit;
RQA.rOuterFit = rOuterFit;
RQA.thickFit = thickFit;
RQA.nYSegments = nYSegments;
RQA.maxSegLen = maxSegLen;
RQA.meanSegLen = meanSegLen;
RQA.smallSegCount = smallSegCount;
RQA.segmentDensity = segmentDensity;
RQA.densityFit = densityFit;
RQA.densityResidual = densityResidual;

%% ===== Angle reference and anatomical sectors =====
thetaDeg = mod(rad2deg(theta),360);
thetaRef = mod(thetaDeg + 60,360);  % adjust anatomical reference

anterior  = thetaRef < 45 | thetaRef >= 315;
rightLat  = thetaRef >= 45  & thetaRef < 135;
posterior = thetaRef >= 135 & thetaRef < 225;
leftLat   = thetaRef >= 225 & thetaRef < 315;

%% ===== Residual analysis =====
outerResidual = rOuterN - rOuterFit;
innerResidual = rInnerN - rInnerFit;
thickResidual = thickN - thickFit;

outerErr = abs(outerResidual);
thickErr = abs(thickResidual);

radial_outer_score = max(0, 100 * (1 - median(outerErr,'omitnan') * 5));
radial_thick_score = max(0, 100 * (1 - median(thickErr,'omitnan') * 5));

%% ===== Segment pattern analysis =====
% posterior zone can naturally have more segments due to psoas/back muscles
multiSeg = nYSegments > 2;

multiSegPosterior = multiSeg & posterior;
multiSegNonPosterior = multiSeg & ~posterior;

posteriorMultiRatio = mean(multiSegPosterior,'omitnan');
nonPosteriorMultiRatio = mean(multiSegNonPosterior,'omitnan');

transition_score = max(0, 100 * (1 - nonPosteriorMultiRatio * 4));

%% ===== Add to output =====
RQA.thetaDeg = thetaDeg;
RQA.thetaRef = thetaRef;

RQA.zone.anterior = anterior;
RQA.zone.rightLat = rightLat;
RQA.zone.posterior = posterior;
RQA.zone.leftLat = leftLat;

RQA.outerResidual = outerResidual;
RQA.innerResidual = innerResidual;
RQA.thickResidual = thickResidual;

RQA.radial_outer_score = radial_outer_score;
RQA.radial_thick_score = radial_thick_score;
RQA.transition_score = transition_score;

RQA.multiSeg = multiSeg;
RQA.multiSegPosterior = multiSegPosterior;
RQA.multiSegNonPosterior = multiSegNonPosterior;
RQA.posteriorMultiRatio = posteriorMultiRatio;
RQA.nonPosteriorMultiRatio = nonPosteriorMultiRatio;

end