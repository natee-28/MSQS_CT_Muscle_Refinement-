function QA = calcMSQS(Image_CT, Mask_Fa, Y, spacing)

Y = logical(Y);
Mask_Fa = logical(Mask_Fa);

% Body wall distance map
distFromSkin = bwdist(~Mask_Fa);

% กำหนด central/deep region ที่ไม่น่าจะเป็น abdominal wall muscle
% ปรับค่า 45-70 pixel ตามภาพ 512x512
central_region = Mask_Fa & distFromSkin > 60;
central_raw = Y & central_region;

% หา component ของ Y
CC = bwconncomp(Y);
L = labelmatrix(CC);

central_bad = false(size(Y));

for i = 1:CC.NumObjects
    comp = (L == i);

    comp_area = sum(comp(:));
    comp_central = sum(comp(:) & central_region(:)) / max(comp_area,1);

    % ถ้า component อยู่ central มาก และมีขนาดเล็ก/แยกเดี่ยว ค่อยถือว่าผิด
    if comp_central > 0.60 && comp_area < 3000
        central_bad = central_bad | comp;
    end
end

wall_band = Mask_Fa & distFromSkin <= 45;

% posterior muscle อาจลึกกว่า ให้เปิดพื้นที่ด้านหลังไว้
[H,W] = size(Image_CT);
posterior_region = false(H,W);
posterior_region(round(H*0.55):end,:) = true;


anatomic_muscle_zone = wall_band | posterior_region;

MuscleCandidate = ...
    (Image_CT > -50) & ...
    (Image_CT < 150) & ...
    Mask_Fa & ...
    anatomic_muscle_zone;

% ถ้า Y เข้า central region มาก ให้ถือว่า leakage
%central_leak = Y & central_region;
central_leak = central_bad;
central_leak_ratio = sum(central_leak(:)) / max(sum(Y(:)),1);
central_leak_score = max(0, 100 * (1 - central_leak_ratio * 10));

MuscleHU = MuscleCandidate; %(Image_CT > -50) & (Image_CT < 150) & Mask_Fa;
organ_leak = Y & ~anatomic_muscle_zone & Mask_Fa;
organ_leak_ratio = sum(organ_leak(:)) / max(sum(Y(:)),1);
organ_leak_score = max(0, 100 * (1 - organ_leak_ratio * 10));

leak_pixels = Y & ~MuscleHU;
leak_ratio = sum(leak_pixels(:)) / max(sum(Y(:)),1);
leak_score = max(0, 100 * (1 - leak_ratio * 5.0));

%Y_dil = imdilate(Y, strel('disk',8));
nearY = imdilate(Y, strel('disk',10)) & ~Y;
raw_missing = MuscleCandidate & nearY;
% ตัดเส้นบาง ๆ รอบ body wall ออก
raw_missing = imopen(raw_missing, strel('disk',2));

% เอาเฉพาะก้อนที่ใหญ่พอ
raw_missing = bwareaopen(raw_missing, 250);

% เอาเฉพาะก้อนที่อยู่ใกล้ Y จริง ๆ
CCm = bwconncomp(raw_missing);
Lm = labelmatrix(CCm);

%missing_pixels = MuscleCandidate & nearY;
missing_pixels = false(size(Y));

for i = 1:CCm.NumObjects
    comp = (Lm == i);

    % ต้องอยู่ติดกับ Y หลัง dilate
    touchY = any((imdilate(comp, strel('disk',3)) & Y), 'all');

    compArea = sum(comp(:));

    if touchY && compArea > 250
        missing_pixels = missing_pixels | comp;
    end
end


%missing_pixels = bwareaopen(missing_pixels,150);

%missing_pixels = MuscleHU & Y_dil & ~Y;
missing_ratio = sum(missing_pixels(:)) / max(sum(Y(:)),1);
%missing_ratio = sum(missing_pixels(:)) / max(sum(MuscleHU & Y_dil,'all'),1);
missing_score = max(0, 100 * (1 - missing_ratio * 1.0));

CC = bwconncomp(Y);
numComp = CC.NumObjects;
shape_score = max(0, 100 - max(0,numComp-10)*5);

hu_vals = Image_CT(Y);
if isempty(hu_vals)
    hu_score = 0;
    meanHU = NaN;
    stdHU = NaN;
else
    meanHU = mean(hu_vals);
    stdHU = std(double(hu_vals));
    outHU_ratio = mean(hu_vals < -50 | hu_vals > 150);
    hu_score = max(0, 100 * (1 - outHU_ratio * 5));
end

perim = bwperim(Y);
areaY = sum(Y(:));
smooth_index = sum(perim(:)) / max(sqrt(areaY),1);
smooth_score = max(0, 100 - max(0, smooth_index - 25)*3);

MSQS = ...
    0.15 * leak_score + ...
    0.15 * missing_score + ...
    0.20 * central_leak_score + ...
    0.20 * organ_leak_score + ...
    0.10 * shape_score + ...
    0.10 * hu_score + ...
    0.10 * smooth_score;

if MSQS >= 85
    action = "ACCEPT";
elseif MSQS >= 65
    action = "REVISE";
else
    action = "MANUAL";
end

QA = struct();
QA.MSQS = MSQS;
QA.action = action;
QA.leak_score = leak_score;
QA.missing_score = missing_score;
QA.shape_score = shape_score;
QA.hu_score = hu_score;
QA.smooth_score = smooth_score;
QA.meanHU = meanHU;
QA.stdHU = stdHU;
QA.leak_pixels = leak_pixels;
QA.missing_pixels = missing_pixels;
QA.leak_ratio = leak_ratio;
QA.missing_ratio = missing_ratio;
QA.num_components = numComp;
QA.central_leak_score = central_leak_score;
QA.central_leak = central_leak;
QA.central_leak_ratio = central_leak_ratio;
QA.organ_leak_score = organ_leak_score;
QA.organ_leak = organ_leak;
QA.organ_leak_ratio = organ_leak_ratio;
QA.MuscleCandidate = MuscleCandidate;

end