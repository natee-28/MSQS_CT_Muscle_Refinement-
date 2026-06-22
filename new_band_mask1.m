NewBandMask = false(size(Y));

for k = 1:numel(RQA.theta)
    t = RQA.theta(k);

    r1 = round(rInnerLocalFull(k) * RQA.bodyR(k));
    r2 = round(RQA.rOuterFit(k)   * RQA.bodyR(k));

    if ~isfinite(r1) || ~isfinite(r2), continue; end
    if r2 < r1, continue; end

    rr = r1:r2;

    x = round(cx + rr .* cos(t));
    y = round(cy + rr .* sin(t));

    valid = x>=1 & x<=size(Y,2) & y>=1 & y<=size(Y,1);
    idx = sub2ind(size(Y), y(valid), x(valid));

    NewBandMask(idx) = true;
end


NewBandMask = imdilate(NewBandMask, strel('disk',1));
NewBandMask = imclose(NewBandMask, strel('disk',2));
NewBandMask = imfill(NewBandMask,'holes');

New_Y = Y & NewBandMask;

%% ........Septum 
% ใช้ NewFit เป็น inner boundary
rSeptumInner = rInnerLocalFull;

% จุดที่ไม่มี muscle / inner เป็น NaN ให้ใช้ outer แทน
bad = ~isfinite(rSeptumInner);
rSeptumInner(bad) = RQA.rOuterFit(bad);

% smooth กันสะดุด
rSeptumInner = movmedian(rSeptumInner,15,'omitnan');

% ถ้ายังมี NaN ใช้ outer เติม
bad = ~isfinite(rSeptumInner);
rSeptumInner(bad) = RQA.rOuterFit(bad);

% safety
rSeptumInner = min(rSeptumInner, RQA.rOuterFit - 0.005);
interpolate_septum_Y = false(size(Y));

for k = 1:numel(RQA.theta)
    t = RQA.theta(k);

    r1 = round(rSeptumInner(k) * RQA.bodyR(k));
    r2 = round(RQA.rOuterFit(k) * RQA.bodyR(k));

    if ~isfinite(r1) || ~isfinite(r2), continue; end
    if r2 < r1, continue; end

    rr = r1:r2;

    x = round(cx + rr .* cos(t));
    y = round(cy + rr .* sin(t));

    valid = x>=1 & x<=size(Y,2) & y>=1 & y<=size(Y,1);
    idx = sub2ind(size(Y), y(valid), x(valid));

    interpolate_septum_Y(idx) = true;
end

interpolate_septum_Y = imdilate(interpolate_septum_Y, strel('disk',1));
interpolate_septum_Y = imclose(interpolate_septum_Y, strel('disk',2));
