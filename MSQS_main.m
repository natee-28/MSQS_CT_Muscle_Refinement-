%%..................Main MSQS  ............
% Date: 4/6/2026...........................................................
% Aimed Plan:
%   1. Open Pts-Folder Dicom .....................................
%   2. Read and sorting DICOM , while kept pixel-data  
%   3. Show as MIP_Coronal or Sagittal to detect L3 or L1(pancreatic-slice)
%   4. Read the slice , then open AI-Model (kept as .mat) 
%   5. Predict Y_pred 
%   6. Post-Process
%   6. MSQS 
%   7. Refine 
%   8. Calculated SMA/SMD and Excel 
%%..............................Natee I.(MI)...............................

%% .............START Here .................................................
clc;clear;
addpath('F:\CT_Sarcopenia_detection\MSQS_Proj')

%%........1. Open Pts-Folder Dicom ................................
uFol = uigetdir();

%%.........

G_col = dicomCollection(uFol) %, "IncludeSubfolders", true, "OutputType","table");
[~,ord] = sort(G_col.StudyDateTime);
G_col = G_col(ord,:);   
%%..4. Read the slice , then open AI-Model (kept as .mat).................
load('model_ct_muscle.mat');
Image_CTal = zeros(512,512,height(G_col));
Meta_al ={};
Yp_all = zeros(size(Image_CTal));
mask_M_al = zeros(size(Image_CTal));
mask_Fa_al = zeros(size(Image_CTal));
Yp_m_all = zeros(size(Image_CTal));


%% .............................................................
for r = 1:1 %height(G_col)
    cycle_ = r;
    files = G_col.Filenames{r};
    [I, meta, spacing] = readSeriesSorted(files);
    fprintf('Volume size: %dx%dx%d  | spacing: [%.3f %.3f %.3f] mm\n', size(I), spacing);
    
    %%..3.Show as MIP_Coronal or Sagittal to detect L3 or L1(pancreatic-slice)
    [n,m,sl] = size(I);
    Mip_I = flipud(squeeze(max(I,[],1)).');  % flipud will change sorting 
    figure(1);imshow(Mip_I,[-500 1700]); axis image; colormap gray
    title('Click the desired slice on MIP');
    
    [x, y] = ginput(1);   %#ok<ASGLU>  % x not used here
    row = round(y);
    
    % Convert clicked row -> slice index (your original idea)
    sL = sl - row;
    
    % safety clamp
    sL = max(1, min(sl, sL));
    close all;
    %%.............wait ......................................................... 
    %sL = sl-cursor_info.Position(2); % Try the last-slice - cursor-slice 
    %........after selected L3 or L1 ........................................ 
    Image_CT = I(:,:,sL);
    imagesc(Image_CT);
    Meta_al{r}=meta(1);
    if (m > n)||(n > m)
    Image_CT = imresize(Image_CT,[512 512],"nearest");
    elseif (n ~= 512);
    Image_CT = imresize(Image_CT,[512 512],"nearest"); 
    end
    %%........Threshold mask .................................................
    Image_CTal(:,:,r)=Image_CT;
    mask_F =   (Image_CT <-30) & (Image_CT > -200);
    Mask_A = imfill((Image_CT > -250) ,'holes'); % For using entire slice with skin 
    
    % Defind All entire slice , and firstly separate M vs Fat.................
    %.mask_fa = Mask fat all...................................................
    CC = bwconncomp(Mask_A);
    F_c1=regionprops(Mask_A,'Area','Centroid');
    idx=find([F_c1.Area]==max([F_c1.Area]));
    center = F_c1(idx).Centroid
    Mask_Fa = ismember(labelmatrix(CC), idx);
    mask_F1 = mask_F .* Mask_Fa;  
    %.. defind Muscle area......................................................
    Mask_aM =  (Image_CT >-50) & (Image_CT < 150);
    mask_M = Image_CT.*Mask_Fa.*Mask_aM;
    
    %%....AI generate Muscle segmentation........................................
    % . 5. Find Abd-Muscle,and using Threshoding to define Muscle and FAT
    % ..test by add envelope around by using edges + image
    % Ypred1 = predict(net,(mask_Me./150))>0.9 ;
    Ypred = predict(net,(mask_M./150))>0.25 ;
    Yp_all(:,:,r) = Ypred;
    % Ha = histogram(nonzeros(mask_M.*Ypred))
    %%...........Muscle defind .............................................  
    figure(1);imagesc(mask_M);hold on;[x,y] = find(Ypred.*(mask_M>0)); plot(y,x,'r.');
    [B_pm,L_pm] = bwboundaries(Ypred);
    figure(2);imshow(squeeze((Image_CT)),[]);
    hold on
    for k = 1:length(B_pm)
    boundary = B_pm{k};
    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2); % waitforbuttonpress;
    end
    close all;

    %cycle_ = [1:r];

%% ...........defind continue..............................................
   %Ypred = logical(Yp_all(:,:,i));
   %mask_M = squeeze(mask_M_al(:,:,i));
   %Mask_Fa = squeeze(mask_Fa_al(:,:,i));
   mean_SMD = mean(nonzeros(mask_M.*Ypred))
   std_SMD = std(nonzeros(mask_M.*Ypred))
   %[Ypred_m] = Improve_pre_M(Image_CT,Mask_Fa,mask_M,net);
   %Yp_m_all(:,:,r) = Ypred;
   %%......check P if not improving using Ypred..............................
   imagesc(Ypred);
   %w = waitforbuttonpress;
   %close all;
     

   Y = bwareaopen(Ypred>0.25, 100);
   QA1 = calcMSQS(Image_CT, Mask_Fa, Y, spacing);
   showMSQSReview(Image_CT, Mask_Fa, Y, QA1);

   if QA1.MSQS >= 85
       Y_final = Y;
       action = "ACCEPT";

   elseif QA1.MSQS >= 65
      [Y_final, refineLog] = refineByActiveContour(Image_CT, Mask_Fa, Y, QA1, spacing);

    else
       Y_final = Y;
       action = "MANUAL";
    end

    QA2 = calcMSQS(Image_CT, Mask_Fa, Y_final, spacing);

    px = spacing(1);
    py = spacing(2);
    SMA = sum(Y_final(:)) * px * py * 0.01;
    SMD = mean(Image_CT(Y_final));
    end
    
    %{
    %% .......semi_automate
    % --- show for QC ---
    figure; imshow(Image_CT.*Y,[-30 150]); axis image; colormap gray
    title('E=ROI refine | G=keep');
    
    waitforbuttonpress;
    key = lower(get(gcf,'CurrentCharacter'));
    
    if key == 'e'
        title('Draw ROI (double-click to finish)');
        cc = roipoly;
        if ~isempty(cc)
            Y = Y & cc;     % refine only
        end
    end
    close(gcf);
    
    
    %% ..check 
    %Ypred = Ypred_m>0.25;
    px = spacing(1);
    py = spacing(2);
    Y_band = (Y).*band;
    [vat_mask,sat_mask] = ray_Gem_casting(center,Y_band,Mask_Fa,Image_CT);
    Ypred = (Y);
    Sub_Fat_area =  sum(logical(sat_mask(:))).*px.*py.*0.01
    Vis_Fat_area = sum(logical(vat_mask(:))).*px.*py.*0.01
    
    %% ........run call arewa .................................................
    %run_cal_area;
    %% ...
    close all;
    p_al = (Ypred.*(mask_M>0));
    n_mus_pixs = sum(p_al(:));
    p_SMA = n_mus_pixs*px.*py*0.01

    %% ... change Sub_Fat_area........................................... 
    imagesc([(20.*out.SAT)+(vat_mask.*25)+(Ypred*5)]);
    Sub_Fat_area = out.SAT_cm2
    %Vis_Fat_area = out.VAT_cm2;
    w = waitforbuttonpress;

    %% ................Sent into Excel .......................................
    % check path and detect test.xls.........................................%%
    close all;
    info = meta(sL);
    cd 'D:\Natee_matlab\Muscle_sacopenia\Jeed_SM_SAT_VAT'; % file xls was here 
    %%...create fields ........................................................%%
    
    a_fields={...
    'PatientID',...
    'gender',...
    'PatientAge',...
    'studydate',...
    'pixel_pm',...
    'SMA',...
    'mean_SMD',...
    'std_SMD',...
    'area_subc_fat',...
    'area_vis_fat',...
    'cycles',...
    'Height',...
    };
    
    %%...............fill calculated parameter ............................%%
    %if  (info(1).PatientAge)
    %    pt_age=info(1).PatientAge;
    %else
    pt_age = str2num(info.StudyDate(1:4))-str2num(info.PatientBirthDate(1:4));
    %end
    % fill the height of pt if u know that 
    prompt = {'Height of Patient: if u know '};
    dlg_title = 'Height of Patient:';
    num_lines = 1;
    def = {'10'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    height_pt = str2num(answer{1});
    close all;
    
    fill_f={info(1).PatientID,...
    info(1).PatientSex,...
    pt_age,...
    info(1).StudyDate,...
    n_mus_pixs,...
    p_SMA,...
    mean_SMD,...
    std_SMD,...
    Sub_Fat_area,...
    Vis_Fat_area,...
    cycle_,...
    height_pt,...
    };  
    %%.........................Fill to Excel ................................%
    if exist('jeed_data_FU.xlsx')
    %check=readtable('test.xls');
    T = readtable('jeed_data_FU.xlsx','ReadRowNames',false);
    check=size(T);
    nn=check(1);
    new_nn=nn+3;
    data_1=table(fill_f);
    formatSpec ='A%d';
    ar=sprintf(formatSpec,new_nn);
    writetable(data_1,'jeed_data_FU.xlsx','Sheet',1,'Range',ar,'WriteVariableNames',false);
    end
    close all; 
    %cd (k_name);
end 
    %}
   


