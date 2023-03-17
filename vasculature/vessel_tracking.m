%% What does this script do?

%% set file path
datapath  = strcat('E:\Jiarui\Data\190426volume\volume\');

% add subfunctions for the script
addpath('E:\Jiarui\Data\code\');

% get the directory of all image tiles
cd(datapath);

tic

% kern=ones(21,21,21)./(21*21*21);

for index=[145:151 168:174 191:197]
    
    v=[];
    % v_low=[];
    % v_high=[];
    
for i=28:31
    load(strcat(num2str(i),'-',num2str(index),'.mat'));
    slice=slice+15;
    slice=depth_corr(slice,2.5);
    % temp=slice(51:150,:,:);
    % temp_low=convn(temp,kern,'same');
    % temp_high=temp-temp_low;
    % temp=temp-min(min(min(temp)));
 
    temp=slice(71:135,:,:);
    % temp_low=temp_low(21:85,:,:);
    % temp_high=temp_high(21:85,:,:);
    v=cat(1,v,temp);
    % v_low=cat(1,v_low,temp_low);
    % v_high=cat(1,v_high,temp_high);
end

%v=v./(max(max(max(v))));
%     for i=1:size(v,1)
%         raw=squeeze(v(i,:,:));
%         tiffname=strcat(num2str(index),'.tif');
%         t = Tiff(tiffname,'a');
%         tagstruct.ImageLength     = 400;%size(data,2);
%         tagstruct.ImageWidth      = 400;%size(data,3);
%         tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
%         tagstruct.BitsPerSample   = 16;
%         tagstruct.SamplesPerPixel = 1;
%         tagstruct.Compression     = Tiff.Compression.None;
%         tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
%         tagstruct.Software        = 'MATLAB';
%         t.setTag(tagstruct);
%         t.write( uint16(65535.*mat2gray(raw)));
%         t.close();
%     end
save([num2str(index),'_depth_corr.mat'],'v');
% save([num2str(index),'_low.mat'],'v_low');
% save([num2str(index),'_high.mat'],'v_high');
toc
end