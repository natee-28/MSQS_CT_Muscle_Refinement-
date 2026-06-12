%{
RQA = radialConsistencyScore(Mask_Fa,Ypred);

figure(1);
plot(rad2deg(RQA.theta), RQA.rOuterN,'.'); hold on;
plot(rad2deg(RQA.theta), RQA.rOuterFit,'LineWidth',2);
xlabel('Angle (degree)');
ylabel('Normalized radius');
title('Outer muscle boundary along radial direction');
residual = RQA.rOuterN - RQA.rOuterFit;
figure(2);
plot(rad2deg(RQA.theta), residual);
%% ..........
[~,idx] = min(residual);
theta_bad = RQA.theta(idx);
S = regionprops(Mask_Fa,'Centroid');
C = S(1).Centroid;
cx = C(1);
cy = C(2);
r = 250;
x2 = cx + r*cos(theta_bad);
y2 = cy + r*sin(theta_bad);
figure;
imshow(Image_CT,[-50 150]);
hold on;
plot([cx x2],[cy y2],'r-','LineWidth',2);
%}
%% ...........................................................................................
RQA = radialConsistencyScore(Mask_Fa,Y_final);
[thetaSort,idx] = sort(RQA.thetaRef);
figure;
plot(thetaSort,RQA.rOuterN(idx),'.'); hold on;
plot(thetaSort,RQA.rOuterFit(idx),'LineWidth',2);
xlabel('Theta ref (degree)');
ylabel('Normalized outer radius');
title('Anatomical-sector radial profile');

figure;
plot(thetaSort,RQA.nYSegments(idx),'.-');
xlabel('Theta ref (degree)');
ylabel('Number of Y segments');
title('Ray segment pattern');

fprintf('Radial outer score: %.1f\n', RQA.radial_outer_score);
fprintf('Thickness score   : %.1f\n', RQA.radial_thick_score);
fprintf('Transition score  : %.1f\n', RQA.transition_score);
fprintf('Posterior multi   : %.3f\n', RQA.posteriorMultiRatio);
fprintf('Non-post multi    : %.3f\n', RQA.nonPosteriorMultiRatio);
fprintf('inner_leak_score   : %.3f\n', RQA.inner_leak_score);
fprintf('inner_miss_score   : %.3f\n', RQA.inner_miss_score);
fprintf('inner_leak_ratio   : %.3f\n', RQA.inner_leak_ratio);
fprintf('inner_miss_ratio   : %.3f\n', RQA.inner_miss_ratio);




%%......................


figure;
plot(thetaSort,RQA.maxSegLen(idx),'.-');
xlabel('Theta ref (degree)');
ylabel('Max segment length');
title('Largest segment length by ray');

figure;
plot(thetaSort,RQA.smallSegCount(idx),'.-');
xlabel('Theta ref (degree)');
ylabel('Small segment count');
title('Small isolated segment pattern');

figure;
plot(thetaSort,RQA.segmentDensity(idx),'.-');
xlabel('Theta ref (degree)');
ylabel('Segment density');
title('Muscle occupancy along ray');
%%.....
%RQA = radialConsistencyScore(Mask_Fa,Y_final);
%[thetaSort,idx] = sort(RQA.thetaRef);

figure;
plot(thetaSort,RQA.thickN(idx),'.-'); hold on;
plot(thetaSort,RQA.densityFit(idx),...
    'LineWidth',2);

xlabel('Theta ref');
ylabel('Thickness Ratio');
title('Muscle Density Signature');
%%...
figure;
plot(thetaSort,...
    RQA.densityResidual(idx),'.-');

yline(0,'k--');
xlabel('Theta ref');
ylabel('Density residual');

%%... 
%[thetaSort,idx] = sort(RQA.thetaRef);

figure;
plot(thetaSort,RQA.innerLeakPx(idx),'.-');
hold on;
yline(-5,'r--');
yline(5,'r--');
xlabel('Theta ref');
ylabel('Inner boundary error (px)');
title('Inner boundary error by angle');
%% ..
RQA.innerLeakIdx = find(RQA.innerLeakRay);
RQA.innerMissIdx = find(RQA.innerMissRay);
figure;
imshow(Image_CT,[-50 150]); axis image; colormap gray; hold on;
visboundaries(Y_final,'Color','g','LineWidth',1);

S = regionprops(Mask_Fa,'Centroid');
C = S(1).Centroid;
cx = C(1); cy = C(2);

% Red = inner leak
for ii = 1:numel(RQA.innerLeakIdx)
    k = RQA.innerLeakIdx(ii);
    t = RQA.theta(k);
    r = RQA.bodyR(k);

    x2 = cx + r*cos(t);
    y2 = cy + r*sin(t);

    plot([cx x2],[cy y2],'r-','LineWidth',1.5);
    text(x2,y2,sprintf('%.0f',RQA.thetaRef(k)),'Color','y');
end

% Blue = inner miss
for ii = 1:numel(RQA.innerMissIdx)
    k = RQA.innerMissIdx(ii);
    t = RQA.theta(k);
    r = RQA.bodyR(k);

    x2 = cx + r*cos(t);
    y2 = cy + r*sin(t);

    plot([cx x2],[cy y2],'b-','LineWidth',1.5);
    text(x2,y2,sprintf('%.0f',RQA.thetaRef(k)),'Color','y');
end

title('Inner anomaly rays using true index');
%% ...
figure;
plot(RQA.thetaDeg,RQA.innerLeakPx,'.-');
hold on;
yline(-5,'r--');
plot(thetaSort,RQA.densityFit(idx).*50,...
    'LineWidth',2);
xlabel('thetaDeg');
ylabel('inner error px');
%% ...
yyaxis left
plot(RQA.thetaDeg,RQA.innerLeakPx,'.-');
ylabel('Inner residual (px)');

yyaxis left
plot(RQA.thetaDeg,RQA.densityFit(idx),'LineWidth',2);
ylabel('Density signature');

densityGrad = abs(gradient(RQA.densityFit));

yyaxis right
plot( RQA.thetaDeg,densityGrad(idx),'g','LineWidth',2)

%% ..
S = regionprops(Mask_Fa,'Centroid');
C = S(1).Centroid;
cx = C(1);
cy = C(2);
figure;
imshow(Image_CT.*Y,[-50 150]); axis image; colormap gray; hold on;
visboundaries(Y_final,'Color','g','LineWidth',1);
susTheta = [23 50];

for tDeg = susTheta
    r = RQA.bodyR(k);
    t = deg2rad(tDeg);

    x2 = cx + r*cos(t);
    y2 = cy + r*sin(t);

    plot([cx x2],[cy y2],...
        'y-','LineWidth',2);

    text(x2,y2,...
        sprintf('%d',tDeg),...
        'Color','y',...
        'FontSize',12,...
        'FontWeight','bold');
end
% Red = inner leak
for ii = 1:numel(RQA.innerLeakIdx)
    k = RQA.innerLeakIdx(ii);
    t = RQA.theta(k);
    r = RQA.bodyR(k);

    x2 = cx + r*cos(t);
    y2 = cy + r*sin(t);

    plot([cx x2],[cy y2],'r-','LineWidth',1.5);
    text(x2,y2,sprintf('%.0f',RQA.thetaRef(k)),'Color','y');
end

%{
RQA.innerLeakIdx = find(RQA.innerLeakRay);
RQA.innerMissIdx = find(RQA.innerMissRay);
S = regionprops(Mask_Fa,'Centroid');
C = S(1).Centroid;
cx = C(1);
cy = C(2);
figure;
imshow(Image_CT.*Y,[-50 150]); axis image; colormap gray; hold on;
visboundaries(Y_final,'Color','g','LineWidth',1);

for ii = 1:numel(RQA.innerLeakIdx)
    k = RQA.innerLeakIdx(ii);
    t = RQA.theta(k);
    r = RQA.bodyR(k);

    x2 = cx + r*cos(t);
    y2 = cy + r*sin(t);

    plot([cx x2],[cy y2],'r-','LineWidth',1.5);
end
for ii = 1:numel(RQA.innerMissIdx)
    k = RQA.innerMissIdx(ii);
    t = RQA.theta(k);
    r = RQA.bodyR(k);

    x2 = cx + r*cos(t);
    y2 = cy + r*sin(t);

    plot([cx x2],[cy y2],'g-','LineWidth',1.5);
end
%}
