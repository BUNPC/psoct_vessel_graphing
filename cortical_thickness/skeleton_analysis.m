%% add path and load file
clear all
addpath('E:\Jiarui\Data\code\');
addpath('E:\Jiarui\Data\code\cortical_thickness');
subject_ID = 'NC_21499';
filepath = ['E:\Jiarui\Data\220401_Ann_15_samples\', subject_ID, '\'];
cd(filepath);
infra = TIFF2MAT('thickness_skel_infra.tif');
supra = TIFF2MAT('thickness_skel_supra.tif');
cortex = TIFF2MAT('thickness_skel_cortex.tif');

%% analyze the thicknes
clc
layer = 'cortex';
roi_sul = ReadImageJROI('RoiSet_sul.zip');
roi_ref = ReadImageJROI('RoiSet_ref.zip');
region_sul = ROIs2Regions(roi_sul,[size(infra,2), size(infra,1)]);
region_ref = ROIs2Regions(roi_ref,[size(infra,2), size(infra,1)]);

thickness_sul = [];
thickness_ref = [];

for i=1:size(cortex,3)
    index_sul = min(ceil(i/5),length(roi_sul));
    index_ref = min(ceil(i/5),length(roi_ref));
    if strcmp(layer,'supra')
        img = transpose(supra(:,:,i));
    elseif strcmp(layer,'infra')
        img = transpose(infra(:,:,i));
    elseif strcmp(layer,'cortex')
        img = transpose(cortex(:,:,i));
    end
    % check if ROI makes sense
%     if i==10
%         img_ref = img;
%         img_sul = img;
%         img_ref(region_sul.PixelIdxList{index_sul})=100;
%         img_sul(region_ref.PixelIdxList{index_ref})=100;
%         figure;imagesc(img_sul);title('sulcus ROI');
%         figure;imagesc(img_ref);title('crest ROI'); 
%     end
    thickness_sul = [thickness_sul; nonzeros(img(region_sul.PixelIdxList{index_sul}))];
    thickness_ref = [thickness_ref; nonzeros(img(region_ref.PixelIdxList{index_ref}))];
end

% calculate statisctics
% mean_sul = mean(thickness_sul,'omitnan')
% std_sul = std(thickness_sul,'omitnan')
% mean_ref = mean(thickness_ref,'omitnan')
% std_ref = std(thickness_ref,'omitnan')
% mean_ref - mean_sul
sqrt((std(thickness_ref,'omitnan')^2+std(thickness_sul,'omitnan')^2)/15)

%% trace the path of skeleton and plot the line profile
% load data
infra_mask = TIFF2MAT([subject_ID '_mask_infra.tif']);
infra_mask = imbinarize(infra_mask);
supra_mask = TIFF2MAT([subject_ID '_mask_supra.tif']);
supra_mask = imbinarize(supra_mask);
cortex_mask = TIFF2MAT([subject_ID '_mask_cortex.tif']);
cortex_mask = imbinarize(cortex_mask);
mus = TIFF2MAT('mus.tif');
layer_class={'infra','supra','cortex'};

% get line profile for each layer
for class = 3%1:length(layer_class)
    layer = layer_class{class};
    i = 14;
    if strcmp(layer,'infra')
        img = squeeze(infra(:,:,(i-1)*5+1));
        mask = squeeze(infra_mask(:,:,i));
    elseif strcmp(layer,'supra')
        img = squeeze(supra(:,:,(i-1)*5+1));
        mask = squeeze(supra_mask(:,:,i));
    elseif strcmp(layer,'cortex')
        img = squeeze(cortex(:,:,(i-1)*5+1));
        mask = squeeze(cortex_mask(:,:,i));
    end
    
    % change mask from infra to supra for different layeres
    mus_img = squeeze(mus(:,:,i)).*single(mask);
    BW = bwskel(imbinarize(img));
    % draw seeding point and the path length wanted
    f1 = figure;
    imagesc(BW);title('Skeleton image (selecting starting point and length)');
    origin = drawpoint;
    path = drawpolyline;
    dx = abs(diff(path.Position(:,1)));
    dy = abs(diff(path.Position(:,2)));
    pathlength = sum(sqrt(dx.^2 + dy.^2));
    x0 = round(origin.Position(1));
    y0 = round(origin.Position(2));
    close(f1);
    % find nearest point in skeleton to this point
    CC = bwconncomp(BW);
    [row, col] = ind2sub(size(BW),CC.PixelIdxList{1});
    [~, idx] = min((row-y0).^2+(col-x0).^2);
    ys = row(idx); xs = col(idx);
    % trace the countour of the skeleton and draw it
    contour = bwtraceboundary(BW,[ys xs],'E',8,round(pathlength),'counterclockwise');
    f1 = figure;
    ax1 = axes('Parent',f1);ax2 = axes('Parent',f1);
    set(ax1,'Visible','off');set(ax2,'Visible','off');
    h1 = imshow(mat2gray(squeeze(mus(:,:,i))),'Colormap',gray,'Parent',ax1);hold on;
    h2 = imshow(mask,'Colormap',gray,'Parent',ax2);
    set(h2, 'AlphaData', 0.3);linkaxes([ax1 ax2],'xy')
    plot(contour(:,2),contour(:,1),'r','LineWidth',2);...
    scatter(x0,y0,20,'g');scatter(contour(end-20,2),contour(end-20,1),20,'y');
    % extract the normal of the skeleton & plot the thickness along the trace of the skeleton
    thickness_profile = [];
    mus_profile = [];
    Orientations = skeletonOrientation(BW,5);
    Normal_ori = Orientations + 90;
    for n = 1:size(contour,1)-20
        thickness_profile(n) = img(contour(n,1),contour(n,2));
        pos_vecs = -round(thickness_profile(n)/2):round(thickness_profile(n)/2);
        angle = Normal_ori(contour(n,1),contour(n,2))/180;
        mus_pos_x = contour(n,2) + cos(angle).*pos_vecs;
        mus_pos_y = contour(n,1) + sin(angle).*pos_vecs;
        mus_pos = sub2ind(size(mus_img), round(mus_pos_y), round(mus_pos_x));
        mus_profile(n) = mean(mus_img(mus_pos));
    end
    % smooth the cortical thickness profile
    thickness_profile = movmean(thickness_profile,5);
    mus_profile = movmean(mus_profile,5);
    cortex_dist = 0.03:0.03:0.03*length(mus_profile);
    figure;yyaxis left;plot(cortex_dist,thickness_profile.*0.03,'LineWidth',1.5);xlabel('Distance from sulcus (mm)');ylabel('cortical thickness (mm)');
        %xlim([0 20]); ylim([0 max(thickness_profile)]);
    yyaxis right;plot(cortex_dist,mus_profile,'LineWidth',1.5);ylabel('scattering coefficient (mm-1)');
        % xlim([0 20]);ylim([0 max(mus_profile)]);title('Measurement along the cortex');
    save(['thickness_', layer, '.mat'], 'thickness_profile');
    save(['mus_', layer, '.mat'], 'mus_profile');
    
end

%% visualization
figure;histogram(thickness_sul,50,'FaceColor','#77AC30');...
    xlabel('Cortical thickness (pxl)');ylabel('# voxels');