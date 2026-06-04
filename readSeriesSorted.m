function [I, meta, spacing] = readSeriesSorted(files)
    meta = arrayfun(@(f) dicominfo(f), files);

    hasIPP = isfield(meta(1),'ImagePositionPatient') && isfield(meta(1),'ImageOrientationPatient');

    if hasIPP
        iop = meta(1).ImageOrientationPatient;
        row = iop(1:3); col = iop(4:6);
        normal = cross(row,col);
        pos = arrayfun(@(m) dot(m.ImagePositionPatient, normal), meta);
        [~,order] = sort(pos);
    else
        if isfield(meta(1),'InstanceNumber')
            inst = arrayfun(@(m) double(m.InstanceNumber), meta);
            [~,order] = sort(inst);
        else
            order = 1:numel(files);
        end
    end

    files = files(order);
    meta  = meta(order);

    nz = numel(files);
    img1 = double(dicomread(files(1)));
    I = zeros(size(img1,1), size(img1,2), nz, "double");
    for k = 1:nz
        I(:,:,k) = double(dicomread(files(k)));
    end

    if isfield(meta(1),'RescaleSlope')
        I = I .* double(meta(1).RescaleSlope) + double(meta(1).RescaleIntercept);
    end

    try
        px = double(meta(1).PixelSpacing(1));
        py = double(meta(1).PixelSpacing(2));
    catch
        px = 1; py = 1;
    end

    dz = NaN;
    if isfield(meta(1),'SpacingBetweenSlices')
        dz = double(meta(1).SpacingBetweenSlices);
    elseif isfield(meta(1),'SliceThickness')
        dz = double(meta(1).SliceThickness);
    elseif hasIPP && numel(meta) > 1
        zLoc = arrayfun(@(m) double(m.ImagePositionPatient(3)), meta);
        dz = median(abs(diff(zLoc)));
    end
    if isnan(dz), dz = 1; end

    spacing = [px py dz];
end