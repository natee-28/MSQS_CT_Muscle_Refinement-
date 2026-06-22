S = regionprops(Mask_Fa,'Centroid');
C = S(1).Centroid;
cx = C(1); cy = C(2);

RQA = radialConsistencyScore(Mask_Fa,Ypred);

cand = find(RQA.nYSegments >= 2);
r2 = [];

for k = cand(:)'
    s = RQA.segStartR{k};
    if numel(s) >= 2
        r2(end+1) = s(2);
    end
end




opts = struct();
opts.rMuscleRefPx = median(r2,'omitnan');
opts.minGapPx = 50;
opts.nearCentroidEndPx = 0.8*opts.rMuscleRefPx;
opts.dilateSeedPx = 2;
opts.extendSingleWindow = 2;
opts.radiusFrac = 0.5;
opts.outerFrac  = 0.90;




OUT = analyzeRQA_Tlevel_Detached(RQA, Y, cx, cy, opts);

figure;
imshow(Image_CT,[-50 150]); axis image; colormap gray; hold on;
visboundaries(Y,'Color','g');
visboundaries(OUT.Yrev,'Color','r','LineWidth',1.5);
title('Detected detached leakage component');

%% ..............
Y_detached_rev=OUT.Yrev; 

RQA = radialConsistencyScore(Mask_Fa,Y_detached_rev);
Y_final = Y_detached_rev; 
debug_plot;
