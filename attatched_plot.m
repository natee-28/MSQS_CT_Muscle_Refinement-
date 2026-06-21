figure;
hold on

plot(RQA.thetaRef,...
     RQA.rInnerFit,...
     'g','LineWidth',2)

plot(RQA.thetaRef,...
     RQA.rOuterFit,...
     'r','LineWidth',2)

plot(RQA.thetaRef,...
     RQA.rOuterFit-RQA.rInnerFit,...
     'k','LineWidth',2)

plot(RQA.thetaRef,...
     RQA.thickN,...
     'b','LineWidth',2)

idx = find(RQA.nYSegments >= 2);

for k = idx(:)'
    xline(RQA.thetaRef(k),'m:');
end

legend('InnerFit',...
       'OuterFit',...
       'ThkFit',...
       'ThkN',...
       'nYS>=2')

xlabel('Theta (deg)')
ylabel('Normalized Radius')
grid on

%% ... 
S = regionprops(Mask_Fa,'Centroid');
C = S(1).Centroid;
cx = C(1);
cy = C(2);
theta = RQA.theta;

xIn = cx + (RQA.rInnerFit .* RQA.bodyR) .* cos(theta);
yIn = cy + (RQA.rInnerFit .* RQA.bodyR) .* sin(theta);

plot(xIn,yIn,'g.','MarkerSize',10)

xOut = cx + (RQA.rOuterFit .* RQA.bodyR) .* cos(theta);
yOut = cy + (RQA.rOuterFit .* RQA.bodyR) .* sin(theta);

plot(xOut,yOut,'r.','MarkerSize',10)

figure;
imshow(Image_CT,[-50 150]);
hold on;

visboundaries(Y,'Color','b');

plot(xIn,yIn,'g.','MarkerSize',10);
plot(xOut,yOut,'r.','MarkerSize',10);
patch([xIn fliplr(xOut)], ...
      [yIn fliplr(yOut)], ...
      'r', ...
      'FaceAlpha',0.4, ...
      'EdgeColor','none');
legend('Y','InnerFit','OuterFit');
%% ...
theta = RQA.theta;

rIn  = RQA.rInnerFit .* RQA.bodyR;
rOut = RQA.rOuterFit .* RQA.bodyR;

xIn  = cx + rIn  .* cos(theta);
yIn  = cy + rIn  .* sin(theta);

xOut = cx + rOut .* cos(theta);
yOut = cy + rOut .* sin(theta);

xPoly = [xIn(:); flipud(xOut(:))];
yPoly = [yIn(:); flipud(yOut(:))];

BandMask = poly2mask(xPoly,yPoly,size(Y,1),size(Y,2));

% optional: intersect with body mask
BandMask = BandMask & Mask_Fa;

Yac = activecontour(Image_CT, BandMask, 1000, 'Chan-Vese');

figure;
imshow(Image_CT,[-50 150]); axis image; colormap gray; hold on;
visboundaries(Y,'Color','g');
visboundaries(BandMask,'Color','y');
visboundaries(Yac,'Color','r');
legend('Y','BandMask','Active contour');