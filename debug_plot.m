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
densityGrad = abs(gradient(RQA.densityFit));
figure;
yyaxis left
plot(RQA.thetaDeg,RQA.innerLeakPx,'.-');
ylabel('Inner residual (px)');

yyaxis right
plot(RQA.thetaDeg,densityGrad(idx),'g','LineWidth',2)
figure;
yyaxis left
plot(RQA.thetaDeg,RQA.densityFit(idx),'LineWidth',2);
ylabel('Density signature');



yyaxis right
plot(RQA.thetaDeg,densityGrad(idx),'g','LineWidth',2)

%% ..
S = regionprops(Mask_Fa,'Centroid');
C = S(1).Centroid;
cx = C(1);
cy = C(2);
figure;
imshow(Image_CT,[-50 150]); axis image; colormap gray; hold on;
visboundaries(Y_final,'Color','g','LineWidth',1);
susidx = [find(RQA.nYSegments>=2)]; 

for ii = 1:numel(susidx)
    k = susidx(ii);
    r = RQA.bodyR(k);
    %t = deg2rad(tDeg);
    t= RQA.theta(k);
   
    x2 = cx + r*cos(t);
    y2 = cy + r*sin(t);

    plot([cx x2],[cy y2],...
        'y-','LineWidth',2);

    text(x2,y2,...
        sprintf('%.0f',RQA.thetaRef(k)),...
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

%% ..
S = regionprops(Mask_Fa,'Centroid');
C = S(1).Centroid;
cx = C(1);
cy = C(2);

figure;
imshow(Image_CT,[-50 150]); axis image; colormap gray; hold on;
visboundaries(Y_final,'Color','g','LineWidth',1);
idx = [find(RQA.nYSegments == 1)];

for ii = 1:numel(idx)
    k = idx(ii);

    t = RQA.theta(k);
    r = RQA.bodyR(k);

    x2 = cx + r*cos(t);
    y2 = cy + r*sin(t);

    plot([cx x2],[cy y2],'r-','LineWidth',1.5);
    text(x2,y2,sprintf('%.0f',RQA.thetaRef(k)),'Color','r');
end
%{
idx = find(~cellfun(@isempty, RQA.segStartR));

for ii = 1:numel(idx)
    k = idx(ii);

    t = RQA.theta(k);
    r = RQA.bodyR(k);

    x2 = cx + r*cos(t);
    y2 = cy + r*sin(t);

    plot([cx x2],[cy y2],'r-','LineWidth',1.5);
    text(x2,y2,sprintf('%.0f',RQA.thetaRef(k)),'Color','r');
end
%}
%% ......
S = regionprops(Mask_Fa,'Centroid');
C = S(1).Centroid;
cx = C(1);
cy = C(2);
figure;
imshow(Image_CT.*Y,[-50 150]); axis image; colormap gray; hold on;
visboundaries(Y_final,'Color','g','LineWidth',1);
idx = [find(RQA.nYSegments >=2)];
for ii = 1:numel(idx)
    k = idx(ii);

    t = RQA.theta(k);
    r = RQA.bodyR(k);

    x2 = cx + r*cos(t);
    y2 = cy + r*sin(t);

    plot([cx x2],[cy y2],'r-','LineWidth',1.5);
    text(x2,y2,sprintf('%.0f',RQA.thetaRef(k)),'Color','r');
end

for k = 1:numel(RQA.segStartR)

    starts = RQA.segStartR{k};
    ends   = RQA.segEndR{k};

    if isempty(starts)
        continue
    end

    t = RQA.theta(k);

    for j = 1:numel(starts)

        r1 = starts(j);
        r2 = ends(j);

        x1 = cx + r1*cos(t);
        y1 = cy + r1*sin(t);

        x2 = cx + r2*cos(t);
        y2 = cy + r2*sin(t);

        plot(x1,y1,'go','MarkerFaceColor','g')
        plot(x2,y2,'bo','MarkerFaceColor','b')

    end

end
%% ... 
S = regionprops(Mask_Fa,'Centroid');
C = S(1).Centroid;
cx = C(1);
cy = C(2);
figure;
%imshow(Image_CT.*Y,[-50 150]); axis image; colormap gray; hold on;
%visboundaries(Y_final,'Color','g','LineWidth',1);
%susTheta = [find(RQA.nYSegments>=2)];

plot(RQA.thetaRef,...
     RQA.rInnerFit.*RQA.bodyR,...
     'k','LineWidth',2)

hold on

for k=1:numel(RQA.segStartR)

    s = RQA.segStartR{k};

    for j=1:numel(s)

        plot(RQA.thetaRef(k),...
             s(j),...
             'ro')
    end

end
%% ..............

Rref = RQA.rInnerFit .* RQA.bodyR;

Dstart = cell(size(RQA.segStartR));

for k = 1:numel(RQA.segStartR)

    s = RQA.segStartR{k};

    if isempty(s) || isnan(Rref(k))
        Dstart{k} = [];
        continue
    end

    Dstart{k} = s - Rref(k);

end

RQA.DstartToBand = Dstart;
idx = find(RQA.nYSegments >= 2);
%idx = find(RQA.nYSegments >= 2);

for ii = 1:numel(idx)
    k = idx(ii);

    segID  = (1:numel(RQA.segStartR{k}))';
    startR = RQA.segStartR{k}(:);
    endR   = RQA.segEndR{k}(:);
    dStart = RQA.DstartToBand{k}(:);

    fprintf('\nRay %d, theta %.1f\n', k, RQA.thetaRef(k));

    T = table(segID, startR, endR, dStart, ...
        'VariableNames', {'Segment','StartR','EndR','Dstart'});

    disp(T)
end
%{
T_remove = 5;  % px, เริ่มต้นก่อน

for k = 1:numel(RQA.segStartR)

    D = RQA.DstartToBand{k};

    if isempty(D)
        continue
    end

    keepSeg   = abs(D) <= T_remove;
    removeSeg = D > T_remove;

end


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
