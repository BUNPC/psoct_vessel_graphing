clear all;

%% ----------------------------------------- %%
% Note Jan 23:

% Current version of code does resampling, background extraction, 
% dispersion compensation, FFT and data stripping, MIP/AIP generation.

% Will implement the surface finding function for data stripping soon.
% Current algorithm needs some further testing.


% - Jiarui Yang
%%%%%%%%%%%%%%%%%%%%%%%

%% set file path
datapath  = strcat('C:\Users\jryang\Downloads\2\');

% add subfunctions for the script
addpath('E:\Jiarui\Data\code\');

%s=zeros(400,400,400);

for iFile=1
% get the directory of all image tiles
    cd(strcat(datapath));
    filename0=dir(strcat('RAW-*.dat'));

    nk = 2048; nxRpt = 2; nx=400; nyRpt = 1; ny = 400;
    dim=[nk nxRpt nx nyRpt ny];
    
    %% image reconstruction
    ifilePath=[datapath,filename0(1).name];
    disp(['Start loading file ', datestr(now,'DD:HH:MM')]);
    [data_ori] = ReadDat_int16(ifilePath, dim); % read raw data: nk_Nx_ny,Nx=nt*nx
    disp(['Raw_Lamda data of file. ', ' Calculating RR ... ',datestr(now,'DD:HH:MM')]);
    data=Dat2RR(data_ori,-0.22);    % tune the second argument to get dispersion compensation
    dat=abs(data(:,:,:));           % calculate the magnitude of the FFT
    
    %% data stripping (crop the depth of the image) and take the log of the image
      slice=dat(1:400,1:400,:);     % will be replaced by the surface finding function
      s=s+slice;
%      mmm=zeros(400);
%      for i=1:400
%          for j=1:400
%              [m, loc]=max(squeeze(slice(21:end,i,j)));
%              mmm(i,j)=loc;
%          end
%      end
%      imagesc(mmm);
%     n=3;
%     v=ones(3,n,n)./(n*n*3);
%     w=convn(slice,v,'same');
%     w=20*log10(w);
    %figure;imagesc(squeeze(slice(:,200,:)));colormap gray;
%     C=strsplit(filename0(iFile).name,'-');
%     coord=strcat(C{3});

    
    % save the average A-line
%      avg=squeeze(mean(w,2));
%      save('avg_white.mat','avg');
%      figure;imagesc(avg);colormap gray;

    %% Generate AIP/MIP and save
    %avg_aline=squeeze(mean(mean(dat,2),3));
    %figure;plot(avg);
    %[m,index]=max(avg_aline);
    % get the tile index

    %aip=squeeze(mean(slice(21:150,:,:),1));
    %mip=squeeze(max(slice(21:end,:,:),[],1));
     %figure;imagesc(aip);colormap gray;axis off;
     %figure;imagesc(squeeze(slice(:,250,:)));title('reconstructed Bscan');ylabel('depth (pxl)');xlabel('X axis (pxl)');colormap gray;
%      figure;imagesc(mip);colormap gray;
%     
%     % define the save path
     %avgname=strcat('C:\Users\jryang\Desktop\Data\0611zf\',num2str(depth),'.mat');
     %mipname=strcat('C:\Users\jryang\Desktop\Data\0611zf\',num2str(depth),'.mat');
     %save(mipname,'mip');
     %save(avgname,'aip');
    
end
        s=s./5;
      matname=strcat(datapath,'data.mat');
      save(matname,'s');
% figure;imagesc(squeeze(slice(140,:,:)));
% index=140;
% img=squeeze(slice(index,:,:));
% img2=squeeze(w(index,:,:));
% 
% m=max(max(img));
% img=uint16(65535*(img./m));
% imwrite(img,'pre5.tif','Compression','none');
% img2=uint16(65535*(img2./m));
% imwrite(img2,'post5.tif','Compression','none');
% avg=mat2gray(avg);
% img=uint16(65535*avg);
% imwrite(img,'bscan_wm.tif','Compression','none');
