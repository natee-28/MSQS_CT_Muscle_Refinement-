function [Y_refined, refineLog] = refineByActiveContour(Image_CT, Mask_Fa, Y, QA, spacing)

Igray = mat2gray(Image_CT, [-50 150]);
Y_refined = logical(Y);

refineLog = struct();
refineLog.mode = "NONE";
refineLog.MSQS_before = QA.MSQS;

figure;
imshow(Image_CT,[-50 150]); axis image; colormap gray; hold on;
visboundaries(Y,'Color','g');
visboundaries(QA.leak_pixels,'Color','r');
visboundaries(QA.missing_pixels,'Color','y');
title('Click missing area to ADD, or press Enter to skip');

[x,y,button] = ginput(1);

if isempty(x)
    refineLog.mode = "SKIP";
    return;
end

cx = round(x);
cy = round(y);

init = false(size(Y));
rad = 20;
r1 = max(1,cy-rad); r2 = min(size(Y,1),cy+rad);
c1 = max(1,cx-rad); c2 = min(size(Y,2),cx+rad);
init(r1:r2,c1:c2) = true;

BW = activecontour(Igray, init, 50, 'Chan-Vese');

MuscleHU = (Image_CT > -50) & (Image_CT < 150) & logical(Mask_Fa);
BW = BW & MuscleHU;

Y_refined = Y_refined | BW;

refineLog.mode = "ADD_ACTIVE_CONTOUR";
refineLog.added_pixels = sum(BW(:));

end