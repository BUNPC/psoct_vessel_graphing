clear all;

%% ----------------------------------------- %%
% Note Jan 23:

% Current version of code does resampling, background extraction, 
% dispersion compensation, FFT and data stripping, MIP/AIP generation.

% Write to TIF images, stitching and blending was done in a seperate script.

% Will implement the surface finding function for data stripping soon.
% Current algorithm needs some further testing.

% Will also integrate write to TIF images, stitching and blending in the
% same script soon.

% - Jiarui Yang
%%%%%%%%%%%%%%%%%%%%%%%

%% set file path
datapath  = strcat('C:\Users\jryang\Google Drive\publication\OCT_fitting\figure\fig4_TDE\depth_profile\');

% add subfunctions for the script
addpath('C:\Users\jryang\Desktop\Data\code\');

% get the directory of all image tiles
cd(datapath);
filename0=dir('1-89*.dat');
mean_I=zeros(400,1);
for iFile=1:length(filename0)
    
    %% add MATLAB functions' path
    %addpath('D:\PROJ - OCT\CODE-BU\Functions') % Path on JTOPTICS
    % addpath('/projectnb/npboctiv/ns/Jianbo/OCT/CODE/BU-SCC/Functions') % Path on SCC server
    
    %% get data information
    nk = 400; nxRpt = 1; nx=400; nyRpt = 1; ny = 400;
    dim=[nk nxRpt nx nyRpt ny];

    %% image reconstruction
    %%%%%%% used for calculating multiple subspectrum, eg. Nspec=[1 2 3 4 8]. set to 1 if calculating only one Nsepc, eg. Nspec=4. 
    ifilePath=[datapath,filename0(iFile).name];
    disp(['Start loading file ', datestr(now,'DD:HH:MM')]);
    [data_ori] = ReadDat_single(ifilePath, dim); % read raw data: nk_Nx_ny,Nx=nt*nx
    data = data_ori(1:400,:,:);
    data = rolloff_corr(data, 2.3);
    mean_I=mean_I+squeeze(mean(mean(data,2),3));
    %aip=squeeze(mean(data_ori(51:250,:,:),1));
%     mip=squeeze(max(slice(201:320,:,:),[],1));
     % figure;imagesc(aip);colormap gray;axis off;
      %figure;plot(squeeze(mean(mean(data,2),3)));ylabel('Intensity (pxl)');xlabel('depth (pxl)');
%     C=strsplit(filename0(iFile).name,'-');
%     coord=strcat(C{3});
%     avgname=strcat('C:\Users\jryang\Desktop\Data\ex_vivo_800um\aip\',coord,'.mat');
% %     mipname=strcat('/projectnb/npbssmic/ns/Jiarui/1219block/mip/vol1/',coord,'.mat');
% %     save(mipname,'mip');
%      save(avgname,'aip');    
%     aip = uint16(65535*(mat2gray(aip)));        % change this line if using mip
%     
%     tiffname=strcat('C:\Users\jryang\Desktop\Data\ex_vivo_800um\aip\',coord,'_aip.tif');
% 
%     t = Tiff(tiffname,'w');
%     tagstruct.ImageLength     = 400;%size(data,2);
%     tagstruct.ImageWidth      = 400;%size(data,3);
%     tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
%     tagstruct.BitsPerSample   = 16;
%     tagstruct.SamplesPerPixel = 1;
%     tagstruct.Compression     = Tiff.Compression.None;
%     tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
%     tagstruct.Software        = 'MATLAB';
%     t.setTag(tagstruct);
%     t.write(aip);
%     t.close();
end
mean_I=mean_I./length(filename0);
% I=mean_I(116:265)-mean(mean_I(end-29:end));
mean_I=mean_I(21:end);
mean_I=10.*log10(mean_I);
figure;plot(mean_I);