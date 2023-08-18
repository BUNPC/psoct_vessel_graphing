dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_10T/';
% Complete subject ID list for Ann_Mckee_samples_10T
subid = {'AD_10382', 'AD_20832', 'AD_20969', 'AD_21354', 'AD_21424',...
         'NC_21499', 'NC_7597',...
         'NC_8095', 'NC_301181'};
subdir = '/dist_corrected/volume/';
for ii =1:length(subid)
    datapath = strcat(dpath, subid{ii}, '/');
    display(datapath)
    %% Generate mask of WM or GM from mus
    vol_mask=[];
    cd(strcat(datapath,'fitting_4x'))
    files=dir('mus*.mat');
    for islice=1:length(files)
        if strcmp(subid{ii},'NC_301181')
            mus_name=strcat('mus',num2str(islice),'_ds4x.mat');
        else
            mus_name=strcat('mus',num2str(islice),'.mat');
        end
        if isfile(mus_name)
            load(mus_name);
            % remove gray matter surface, which has abnormally high mus due to
            % unflat surface
            mask2=zeros(size(MosaicFinal));
            mask2(MosaicFinal>0.01)=1;
            mask2=imerode(mask2,strel('disk',10));
            % dynamically find mus threshold that separates white-gray matter
            groups=0:0.5:25;
            V=histcounts(MosaicFinal(mask2==1),groups);
            L=islocalmax(V);
            l=islocalmin(V);
            WM_center_mus=groups(find(L,1,'last'));
            WG_bound_mus=groups(find(l,1,'last'));
            mask=zeros(size(MosaicFinal));
            % apply threshold to get only WM or GM
%             mask(MosaicFinal>(WM_center_mus+WG_bound_mus)/2)=1; % WM
            mask(MosaicFinal<(WM_center_mus+WG_bound_mus)/2)=1; % GM
            % clean the mask, remove small ilands, fill in small holes
            mask=mask.*mask2;
            mask=imerode(mask,strel('disk',2));
            % dilate the mask to remove empty holes inside tissue, such as
            % in vessels
            mask=imdilate(mask,strel('disk',10));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % remove unconnected ilands, leave only ilands with area larger than 500 pixels
            mask=KeepMajorMask(mask,500);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % generate volume mask for 150um thickness
            for kk =islice*11-10:islice*11
                vol_mask(:,:,kk)=mask;
            end
        end
    end
        
    %% load graph
    total_vessel_length=0;
    cd(strcat(datapath,'dist_corrected/volume'))
    load('ref_4ds_norm_inv_segment_sigma1_thresh0.24_mask40_graph.mat')
    nodes=Graph.nodes;
    edges=Graph.edges;
    % loop over all edges
    for iedge =1:length(edges)
        edge=edges(iedge,:);
        node1=nodes(edge(1),:);
        node2=nodes(edge(2),:);
        % calculate center mass of this edge
        center_mass=round((node1+node2)/2);
        % WM mask from mus is usually thinner in depth than volume recon
        if center_mass(3)<size(vol_mask,3)
            % add length of this edge if its within the WM or GM mask
            if vol_mask(center_mass(2),center_mass(1),center_mass(3))==1
                total_vessel_length=total_vessel_length+sqrt(sum(((node1-node2).*[0.012, 0.012, 0.015]).^2));
            end
        end
    end
    % calculate total vessel length over WM or GM volume
    GM_density(ii)=total_vessel_length/(sum(vol_mask(:)).*0.012*0.012*0.015)
    
end
figure;
bar(1:9,GM_density,'FaceColor',[0,0,0]);
ylabel('vessel total length over GM volume')
title('Gray matter')
ax=gca;
ax.FontSize=15;
xticklabels(subid)
xtickangle(60)
