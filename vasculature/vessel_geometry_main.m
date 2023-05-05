%% Perform geometric analysis of pruned graph.
% This script should be run after running vessel segmentation, converting
% the segment to graph, and pruning the graph with vesGraphValidate.
%{
This script performs the following:
- topology (why does this use the mask?)
- diameter
- length
- tortuosity (vessel_tortuosity_index.m)
%}

%% extract vessels and clean up boundaries by creating a mask
V_seg=TIFF2MAT('I_seg_Ann.tif');
mask=TIFF2MAT('I_mask.tif');
for i=1:size(img,3)
    mask_tmp=squeeze(mask(:,:,i));
    mask_tmp(mask_tmp(:)~=0)=1;
    mask_tmp = bwareaopen(logical(mask_tmp), 500);
    V_seg(:,:,i)=uint16(mask_tmp).*squeeze(V_seg(:,:,i));
    % figure;imshow(mask_tmp);
end
MAT2TIFF(V_seg,'I_seg_masked2.tif');

%% calculate topology
% skeletonization
I_seg=TIFF2MAT('I_seg_nor_mask.tif');
I_seg(I_seg(:)~=0)=1;
I_skel=bwskel(logical(I_seg));
MAT2TIFF(I_skel,'I_seg_skel.tif');

%% calculate radius
global Data;
Diam=GetDiam_graph(permute(Data.angio,[2 3 1]),Data.Graph.nodes,Data.Graph.edges,0.5,[10 10 10]);
dia_vessel=zeros(1,length(Data.Graph.segInfo.segLen));
for i=1:length(dia_vessel)
    dia_vessel(i)=median(Diam(find(Data.Graph.segInfo.nodeSegN(:)==i)));
end

figure;histogram(dia_vessel,'BinWidth',10);


%% length of each vessel
length_vessel=sum(:,1).*sum(:,6);
length_vessel(length_vessel(:)==0)=[];    % get rid of fake vessels
figure;histogram(length_vessel.*10,0:25:1000);

%% Calculate Diameter
%{
% Load Graph struct
fname = 'volume_nor_inverted_masked_sigma1_frangi_seg_regraphed';
dpath =...
    'C:\Users\mack\Documents\BU\Boas_Lab\psoct_human_brain_resources\test_data\Hui_Frangi_dataset\200726PSOCT\';
fpath = strcat(dpath, strcat(fname,'.mat'));
Data = load(fpath,'Graph');

% Call function to calculate diameter at each node
Diam = GetDiam_graph(...
    vol,...
    Data.Graph.nodes,...
    Data.Graph.edges,...
    Ithresh,...
    vox_dim);
%}
%% Calculate Tortuosity
%{
tortuosity = vessel_tortuosity_index(Data.Graph, Ithresh);
%}
%% Histograms for geometries
%{
% Histo for diameter
figure; histogram(Diam);
title('Vessel Diameter')
xlabel('Diameter (microns)')
ylabel('Count')
set(gca, 'FontSize', 20)

% Histo for diameter
figure; histogram(tortuosity);
title('Vessel Tortuosity')
xlabel('Tortuosity (unitless)')
ylabel('Count'); ylim([0,200])
set(gca, 'FontSize', 20)
%}

%% (OLD CODE) Median diameter
%{
dia_vessel=zeros(1,length(Data.Graph.segInfo.segLen));
for i=1:length(dia_vessel)
    dia_vessel(i)=median(Diam(find(Data.Graph.segInfo.nodeSegN(:)==i)));
end
figure; histogram(dia_vessel,'BinWidth',10);
%}