thetaRange = 105:144;

leftRef  = 137:144;
midLeak  = 116:136;
rightRef = 105:115;

rSample = nan(size(RQA.thetaRef));

% 1) ใช้ rInner เดิมจาก nYS==1 ด้านซ้าย
for k = leftRef
    if RQA.nYSegments(k)==1
        rSample(k) = RQA.rInnerN(k);
    end
end
% 2) ช่วง leak ใช้ segment2 เป็น presumed true muscle
for k = midLeak
    if RQA.nYSegments(k) >= 2
        s = RQA.segStartR{k};

        if numel(s) >= 2
            rSample(k) = s(2) ./ RQA.bodyR(k);
        end
    end
end


% 3) ใช้ rInner เดิมจาก nYS==1 ด้านขวา
for k = rightRef
    if RQA.nYSegments(k)==1
        rSample(k) = RQA.rInnerN(k);
    end
end

% fit ใหม่เฉพาะช่วง

x = thetaRange;
y = rSample(thetaRange);

rInnerLocal = y;
rInnerLocal = fillmissing(rInnerLocal,'pchip');
rInnerLocal = movmedian(rInnerLocal,15,'omitnan');

rInnerLocalOnly = nan(size(RQA.thetaRef));
rInnerLocalOnly(thetaRange) = rInnerLocal;

rInnerLocalFull = RQA.rInnerN;
%  Fill  rInnerN only nonfinite section with rInnerFit
bad = ~isfinite(rInnerLocalFull);
rInnerLocalFull(bad) = RQA.rInnerFit(bad);

rInnerLocalFull(thetaRange) = rInnerLocal;
% safety: inner ต้องอยู่ด้านใน outer เสมอ
epsBand = 0.01;
rInnerLocalFull = min(rInnerLocalFull, RQA.rOuterFit - epsBand);

figure; hold on
plot(RQA.thetaRef, RQA.rInnerFit,'g--','LineWidth',1.5)
plot(RQA.thetaRef, rInnerLocalOnly,'r','LineWidth',2)
plot(thetaRange, rSample(thetaRange),'bo','MarkerFaceColor','b')
legend('Original rInnerFit','Local reconstructed','Samples')
grid on

%% ....test again after revised ........................... 
figure;
hold on
plot(RQA.thetaRef,...
     rInnerLocalFull,...
     'k','LineWidth',2)

plot(RQA.thetaRef,...
     RQA.rInnerFit,...
     'g','LineWidth',2)

plot(RQA.thetaRef,...
     RQA.rInnerN,...
     'm','LineWidth',2)

plot(RQA.thetaRef,...
     RQA.rOuterFit,...
     'r','LineWidth',2)

plot(RQA.thetaRef,...
     RQA.rOuterFit-RQA.rInnerFit,...
     'k','LineWidth',2)

plot(RQA.thetaRef,...
     RQA.rOuterFit-rInnerLocalFull,...
     'm','LineWidth',2)

plot(RQA.thetaRef,...
     RQA.thickN,...
     'b','LineWidth',2)

idx = find(RQA.nYSegments >= 2);

for k = idx(:)'
    xline(RQA.thetaRef(k),'m:');
end

legend('newFit',...
       'InnerFit',...
       'InnerN',...
       'OuterFit',...
       'ThkFit',...
       'newThkN',...
       'ThkN',...
       'nYS>=2')

xlabel('Theta (deg)')
ylabel('Normalized Radius')
grid on

%% .......Anatomical map check 

S = regionprops(Mask_Fa,'Centroid');
C = S(1).Centroid;
cx = C(1);
cy = C(2);

theta = RQA.theta;

xIn = cx + (rInnerLocalFull .* RQA.bodyR) .* cos(theta);
yIn = cy + (rInnerLocalFull .* RQA.bodyR) .* sin(theta);

plot(xIn,yIn,'g.','MarkerSize',10)

xOut = cx + (RQA.rOuterFit .* RQA.bodyR) .* cos(theta);
yOut = cy + (RQA.rOuterFit .* RQA.bodyR) .* sin(theta);

plot(xOut,yOut,'r.','MarkerSize',10)
figure; 
plot(xIn,yIn,'g.','MarkerSize',10)
figure;
imshow(Image_CT,[-50 150]);
hold on;

visboundaries(Y,'Color','b');
% .........Gen new 
% ensure row vectors
xIn  = xIn(:)';
yIn  = yIn(:)';
xOut = xOut(:)';
yOut = yOut(:)';

% close loop explicitly
xPoly = [xIn, fliplr(xOut), xIn(1)];
yPoly = [yIn, fliplr(yOut), yIn(1)];

valid = isfinite(xPoly) & isfinite(yPoly);
xPoly = xPoly(valid);
yPoly = yPoly(valid);

NewBandMask = poly2mask(xPoly, yPoly, size(Y,1), size(Y,2));
New_Y = Y & NewBandMask;
%% ...Gen New_Y ...........................................................
xPoly = [xIn(:); flipud(xOut(:))];
yPoly = [yIn(:); flipud(yOut(:))];

valid = isfinite(xPoly) & isfinite(yPoly);
xPoly = xPoly(valid);
yPoly = yPoly(valid);

NewBandMask = poly2mask(xPoly, yPoly, size(Y,1), size(Y,2));
New_Y = Y & NewBandMask;

plot(xIn,yIn,'g.','MarkerSize',10);
plot(xOut,yOut,'r.','MarkerSize',10);
patch([xIn fliplr(xOut)], ...
      [yIn fliplr(yOut)], ...
      'r', ...
      'FaceAlpha',0.4, ...
      'EdgeColor','none');
legend('NewFit','outFit','M_Fit');
