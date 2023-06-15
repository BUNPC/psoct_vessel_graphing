function grayscale = FindGrayScale(file,modality)
%% set path
if isa(file, 'numeric')
    I0 = file;
elseif isa(file,"string") || isa(file,"char")
    load(file,'MosaicFinal');
    I0 = MosaicFinal;
end

% modality = Mosaic2D.Modals;
% path = Mosaic2D.outdir;
% 
% load([path filesep modality{ii} '_slice001.mat']);

%% Load image and find histogram peak
img_e = I0(I0~=0);

nbin = 100;
h = histogram(img_e,nbin);
L = islocalmax(h.Values);
peak = h.BinEdges(L);
range1_99 = [prctile(img_e,1), prctile(img_e,99)];
peak(peak<range1_99(1)|peak>range1_99(2))=[];
xline(peak);

peak_normalized = mat2gray(double(peak),double(range1_99));
npeak = sum(diff(peak_normalized) > 0.07)+1;

fprintf('%i peak/peaks detected.\n',npeak);


%% Generating GrayScale Option
if npeak==1                             % likely from TDE experiment
    switch modality
        case {'AIP','MIP'}
            thre(1) = prctile(img_e,2.5);
            high = prctile(img_e,99);
            Option1 = [thre(1) high];
        case 'Retardance'
            low = prctile(img_e,1);
            thre(1) = prctile(img_e,97.5);
            Option1 = [low thre(1)];
    end
    
    
    
elseif npeak>1
    switch modality
        case {'AIP','MIP'}
            th      = multithresh(img_e,2);
            thre(1) = th(1);             % agar threshold
            high    = prctile(img_e,99.75); % overexposure threshold
            Option1 = [thre(1) high];

            th      = multithresh(img_e,1);
            thre(2) = th(1);             % agar threshold
            high    = prctile(img_e,99.75); % overexposure threshold
            Option2 = [thre(2) high];

        
        case 'Retardance'
            th      = multithresh(img_e,2);
            low     = prctile(img_e,1); 
            thre(1) = th(2);             % agar threshold
            Option1 = [low thre(1)];
            
            th      = multithresh(img_e,1);
            low     = prctile(img_e,1); 
            thre(2) = th(1);             % agar threshold
            Option2 = [low thre(2)];
    end
end

figure; 
subplot(121); histogram(img_e);xline(Option1); 
title(sprintf('Option1 grayscale = [%f %f]',Option1(1),Option1(2)));
subplot(222); imshow(I0, Option1);
subplot(224); imshow(I0<thre(1));
pause(0.1)

if exist('Option2','var')
    figure; 
    subplot(121); histogram(img_e);xline(Option2); 
    title(sprintf('Option2 grayscale = [%f %f]',Option2(1),Option2(2)));
    subplot(222); imshow(I0, Option2);
    subplot(224); imshow(I0<thre(2));
end


%% grayscale selection
answer1 = input('Option1 or Option2 or Customize? [ 1 / 2 / c ]:   ','s');
switch answer1
    case '1'; grayscale = Option1;
    case '2'; grayscale = Option2;
    case 'c'; chosen = false;
end

%% grayscale customization

while chosen == false
    grayscale = input('Please input grayscale value. e.g. [55, 90]:   ');
    figure;
    subplot(121); histogram(img_e);xline(grayscale); title(sprintf('grayscale = [%f %f]',grayscale(1),grayscale(2)));
    subplot(222); imshow(I0, grayscale);
    switch modality
        case {'AIP','MIP'}; subplot(224); imshow(I0<grayscale(1));
        case 'Retardance'; subplot(224); imshow(I0<grayscale(2));
    end
    
    answer2 = input('save result? [ Y / N ]:   ','s'); 
    switch answer2
        case 'Y'; chosen = true;
        case 'N'; chosen = false;
    end
end
end
% air - agar - wm - gm
% tissue vs no tissue
% TDE vs no TDE