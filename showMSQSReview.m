function showMSQSReview(Image_CT, Mask_Fa, Y, QA)

figure;
imshow(Image_CT,[-50 150]); axis image; colormap gray; hold on;

visboundaries(Y,'Color','g','LineWidth',1);
visboundaries(QA.leak_pixels,'Color','r','LineWidth',1);
visboundaries(QA.missing_pixels,'Color','y','LineWidth',1);
visboundaries(QA.organ_leak,'Color','m','LineWidth',1);
title(sprintf('MSQS %.1f | %s | G=AI, R=HU Leak, Y=Missing, C=Central, M=Organ', ...
    QA.MSQS, QA.action));

end