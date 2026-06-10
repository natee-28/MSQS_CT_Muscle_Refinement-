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
%% ......................
RQA = radialConsistencyScore(Mask_Fa,Y_final);
[thetaSort,idx] = sort(RQA.thetaRef);

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
%% .....
RQA = radialConsistencyScore(Mask_Fa,Y_final);
[thetaSort,idx] = sort(RQA.thetaRef);

figure;
plot(thetaSort,RQA.thickN(idx),'.-'); hold on;
plot(thetaSort,RQA.densityFit(idx),...
    'LineWidth',2);

xlabel('Theta ref');
ylabel('Thickness Ratio');
title('Muscle Density Signature');
%% ...
figure;
plot(thetaSort,...
    RQA.densityResidual(idx),'.-');

yline(0,'k--');
xlabel('Theta ref');
ylabel('Density residual');